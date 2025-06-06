/*
 * Copyright (c) 2018. Real Time Genomics Limited.
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

package com.rtg.variant.bayes.complex;

import com.rtg.launcher.globals.CoreGlobalFlags;
import com.rtg.launcher.globals.GlobalFlags;
import com.rtg.mode.DNA;
import com.rtg.util.MaxShiftUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.machine.MachineType;
import com.rtg.variant.MachineErrorChooserInterface;
import com.rtg.variant.ThreadingEnvironment;
import com.rtg.variant.VariantAlignmentRecord;
import com.rtg.variant.VariantParams;
import com.rtg.variant.bayes.Evidence;
import com.rtg.variant.bayes.Factor;
import com.rtg.variant.bayes.Hypotheses;
import com.rtg.variant.bayes.snp.HypothesesPrior;
import com.rtg.variant.match.AlignmentMatch;
import com.rtg.variant.realign.AlignmentEnvironment;
import com.rtg.variant.realign.AlignmentEnvironmentCG;
import com.rtg.variant.realign.AlignmentEnvironmentGenomeSubstitution;
import com.rtg.variant.realign.AlignmentEnvironmentRead;
import com.rtg.variant.realign.AllPaths;
import com.rtg.variant.realign.Environment;
import com.rtg.variant.realign.EnvironmentCombined;
import com.rtg.variant.realign.InvertCgTemplateEnvironment;
import com.rtg.variant.realign.RealignParams;
import com.rtg.variant.realign.ScoreFastUnderflow;
import com.rtg.variant.realign.ScoreFastUnderflowCG;
import com.rtg.variant.util.arithmetic.PossibilityArithmetic;

/**
 * Provides evidence for a complex hypothesis by performing an all-paths alignment
 * from an alignment record to each hypothesis.
 */
public final class EvidenceComplex extends Evidence {

  /** Print complex evidence scores into the developer log for debugging. */
  private static final boolean PRINT_EVIDENCE_DETAILS = GlobalFlags.isSet(CoreGlobalFlags.COMPLEX_EVIDENCE_DETAILS);

  private static final ScoreInterfaceMemoInterface SCORE_INTERFACE_MEMO;

  // If true, CG allpaths realignment should use the full reconstructed read, otherwise use the flattened representation
  static final boolean CG_ALLPATHS = GlobalFlags.getBooleanValue(CoreGlobalFlags.COMPLEX_CALLER_UNROLL_CG_FLAG);

  static {
    SCORE_INTERFACE_MEMO = new ScoreInterfaceMemo();
  }

  private static AllPaths getAllPaths(final RealignParams params) {
    if (CG_ALLPATHS && params.machineType() != null && params.machineType().isCG()) {
      return new ScoreFastUnderflowCG(params);
    }
    return new ScoreFastUnderflow(params);
  }

  private final int mReference;

  private final double[] mProb;

  private double mError = 0.0;

  private double mPE;

  private final Factor<DescriptionComplex> mHypotheses;

  private final int mReadHypothesis;

  private final int mReadBasesLeftOfMatch;
  private final int mReadBasesRightOfMatch;

  private final PossibilityArithmetic mArithmetic;

  private final double mLogSum;

  private static volatile boolean sIsIonTorrent = false;

  private final AlignmentMatch mMatch;

  // This is not very nice
  private static synchronized void printIonTorrentDevLog() {
    if (!sIsIonTorrent) {
      Diagnostic.developerLog("Ion Torrent machine type - widening allpaths");
      sIsIonTorrent = true;
    }
  }

  private final boolean mForward;
  private final boolean mReadPaired;
  private final boolean mFirst;
  private final boolean mMated;
  private final boolean mUnmapped;

  /**
   * @param hypotheses description of the underlying hypotheses.
   * @param match Match object to get distribution for
   * @param reference complex reference
   * @param params variant params
   * @param chooser machine error chooser
   */
  public EvidenceComplex(HypothesesPrior<DescriptionComplex> hypotheses, AlignmentMatch match, ComplexTemplate reference, VariantParams params, MachineErrorChooserInterface chooser) {
    super(hypotheses.description(), match.mapError());
    mHypotheses = hypotheses;
    mReference = hypotheses.reference();
    mArithmetic = hypotheses.arithmetic();
    mPE = mArithmetic.zero();
    mMatch = match;
    mReadBasesLeftOfMatch = match.getBasesLeftOfMatch();
    mReadBasesRightOfMatch = match.getBasesRightOfMatch();
    final int size = description().size();
    mProb = new double[size];
    final VariantAlignmentRecord alignmentRecord = match.alignmentRecord();
    final RealignParams me = chooser.realignParams(alignmentRecord.getReadGroup(), alignmentRecord.isReadPaired());
    final boolean cg = me.machineType() != null && me.machineType().isCG();
    final AllPaths sm;
    if (params.threadingEnvironment() == ThreadingEnvironment.PARALLEL) {
      sm = getAllPaths(me);
    } else {
      sm = SCORE_INTERFACE_MEMO.getScoreInterface(me);
    }

    final AlignmentEnvironment se;
    if (cg) {
      se = new AlignmentEnvironmentCG(alignmentRecord, params, reference.templateBytes(), me.machineType());
    } else {
      if (params.complexUseSoftClip()) {
        se = new AlignmentEnvironmentRead(alignmentRecord, params, me.machineType());
      } else {
        se = new AlignmentEnvironmentRead(alignmentRecord, params, me.machineType(), match.getSoftClipLeft(), match.getSoftClipRight());
      }
    }
    final int maxShift0;
    final int newStart;
    int softClipStartOffset = 0;
    if (me.machineType() == MachineType.IONTORRENT) {
      final MaxShiftCigarParser calc = new MaxShiftCigarParser();
      calc.parse(alignmentRecord.getCigar(), alignmentRecord.getStart());
      maxShift0 = calc.getMaxShift();
      newStart = calc.getStartPos();
      if (params.complexUseSoftClip()) {
        softClipStartOffset = calc.getSoftClipStartOffset();
      }
      if (!sIsIonTorrent) {
        printIonTorrentDevLog();
      }
    } else {
      maxShift0 = MaxShiftUtils.calculateDefaultMaxShift(se.subsequenceLength());
      newStart = alignmentRecord.getStart();
      if (params.complexUseSoftClip()) {
        softClipStartOffset = match.getSoftClipLeft();
      }
    }
    //if (softClipStartOffset > 0) {
    //  System.err.println("EvidenceComplex with soft clipped start: " + alignmentRecord.getStart() + " " + alignmentRecord.getEnd() + " " + alignmentRecord.getCigar());
    //}
    final int adjust = hypotheses.description().maxLength() - hypotheses.description().minLength();
    assert adjust >= 0;

    final int maxShift = maxShift0 + adjust;
    final double[] logScore = new double[size];
    double sum = mArithmetic.zero();
    for (int i = 0; i < size; ++i) {
      //TODO put in fast delta scoring for non-CG case.
      final String replace = description().name(i);
      final AlignmentEnvironment temEnv = new AlignmentEnvironmentGenomeSubstitution(se.start() - softClipStartOffset, 0 /* doesn't matter */, reference, DNA.stringDNAtoByte(replace));
      final EnvironmentCombined envTmp = new EnvironmentCombined(se, newStart - softClipStartOffset, maxShift, temEnv);
      final Environment env;
      if (cg && se.isInverted()) {
        env = new InvertCgTemplateEnvironment(envTmp, me.machineType());
      } else {
        env = envTmp;
      }
      sm.setEnv(env);
      final double poss = mArithmetic.ln2Poss(sm.totalScoreLn());
      //System.err.println("Read match=" + match.readString() + " Hyp i=" + i + " name=" + hypotheses.description().name(i) + " : unnorm score=" + mArithmetic.poss2Ln(poss) + " scorematrix=\n" + sm.toString());
      logScore[i] = poss;
      sum = mArithmetic.add(sum, poss);
    }

    // Normalize, and determine readHyp
    mLogSum = mArithmetic.poss2Ln(sum);
    int readHyp = Hypotheses.NO_HYPOTHESIS;
    double maxProb = -1;
    for (int i = 0; i < size; ++i) {
      final double poss = mArithmetic.divide(logScore[i], sum);
      set(i, poss);
      if (mProb[i] > maxProb) {
        readHyp = i;
        maxProb = mProb[i];
      }
    }
    if (PRINT_EVIDENCE_DETAILS) {
      for (int i = 0; i < size; ++i) {
        final double poss = mArithmetic.divide(logScore[i], sum);
        Diagnostic.developerLog("CX_MATCH: " + (match.isFixedLeft() ? "" : "~") + match.readString() + (match.isFixedRight() ? "" : "~")
          + " hyp: " + (i == readHyp ? "*" : " ") + (i == mReference ? "= " : "X ") + i + " " + hypotheses.description().name(i)
          + " score: " + mArithmetic.poss2Prob(poss)
          + (match.alignmentRecord().getReadGroup() == null ? "" : " sample: " + match.alignmentRecord().getReadGroup().getSample())
          + " cigar: " + match.alignmentRecord().getCigar());
      }
    }

    if (readHyp == Hypotheses.NO_HYPOTHESIS) {
      mReadHypothesis = readHyp;
      mError = 1.0;
    } else {
      mReadHypothesis = mProb[readHyp] > 0.5 ? readHyp : Hypotheses.NO_HYPOTHESIS;
      mError = 1.0 - maxProb;
    }

    mForward = !match.alignmentRecord().isNegativeStrand();
    mReadPaired = match.alignmentRecord().isReadPaired(); // True if is a paired end read at all, false if single end
    mFirst = match.alignmentRecord().isFirst();
    mMated = match.alignmentRecord().isMated(); // True if is a paired end read and is properly mated, false if single end or unmated
    mUnmapped = match.alignmentRecord().isUnmapped();
  }

  AlignmentMatch match() {
    return mMatch;
  }

  /**
   * Get the raw unnormalized sum in log probability.
   * Handy for unit tests.
   * @return the log sum.
   */
  double sumLn() {
    return mLogSum;
  }

  /**
   * Set the distribution value for a specified hypothesis.
   * @param index selects the (haploid) hypothesis.
   * @param poss probability to be assigned (in possibility arithmetic).
   */
  private void set(final int index, final double poss) {
    assert mArithmetic.isValidPoss(poss);
    mProb[index] = mArithmetic.poss2Prob(poss);
    final double hypPrior = mHypotheses.p(index);
    assert mArithmetic.isValidPoss(hypPrior);
    mPE = mArithmetic.add(mPE, mArithmetic.multiply(hypPrior, poss));
    assert mArithmetic.isValidPoss(mPE);
  }

  @Override
  public double probability(int index) {
    return mProb[index];
  }

  @Override
  public double pe() {
    return mArithmetic.poss2Prob(mPE);
  }

  @Override
  public int read() {
    return mReadHypothesis;
  }

  @Override
  public int getReadBasesLeft() {
    return mReadBasesLeftOfMatch;
  }

  @Override
  public int getReadBasesRight() {
    return mReadBasesRightOfMatch;
  }

  @Override
  public void setReadBasesLeft(int readBasesLeft) {
    throw new UnsupportedOperationException();
  }

  @Override
  public void setReadBasesRight(int readBaseRight) {
    throw new UnsupportedOperationException();
  }

  @Override
  public double error() {
    return mError;
  }

  @Override
  public boolean isForward() {
    return mForward;
  }

  @Override
  public boolean isReadPaired() {
    return mReadPaired;
  }

  @Override
  public boolean isFirst() {
    return mFirst;
  }

  @Override
  public boolean isMated() {
    return mMated;
  }

  @Override
  public boolean isUnmapped() {
    return mUnmapped;
  }
}
