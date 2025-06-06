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
package com.rtg.alignment;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Locale;

import com.rtg.launcher.globals.GlobalFlags;
import com.rtg.launcher.globals.CoreGlobalFlags;
import com.rtg.mode.ProteinScoringMatrix;
import com.rtg.ngs.NgsOutputParams;
import com.rtg.ngs.NgsParams;
import com.rtg.protein.GotohProteinEditDistance;
import com.rtg.reader.CgUtils;
import com.rtg.reader.PrereadType;
import com.rtg.reader.SequencesReader;
import com.rtg.util.InvalidParamsException;
import com.rtg.util.MaxShiftFactor;
import com.rtg.util.StringUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.machine.MachineType;
import com.rtg.variant.MachineErrorParamsBuilder;
import com.rtg.variant.realign.RealignParamsImplementation;

import htsjdk.samtools.SAMReadGroupRecord;

/**
 * Factory which creates a specific Edit Distance implementation depending on
 * the parameters passed.
 *
 */
public final class EditDistanceFactory {

  //  private static final String LB_VALUE = "lower-bound-value";

  // Only use the gotoh aligner (disable all others, has priority over USE_NOINDELS_ONLY)
  private static final boolean USE_GOTOH_ONLY = GlobalFlags.isSet(CoreGlobalFlags.EDIT_DIST_GOTOH_ONLY_FLAG);
  // Only use the single-indel-seeded aligner (disable all others, has priority over USE_NOINDELS_ONLY)
  private static final boolean USE_SINGLE_INDEL_SEEDED_ONLY = GlobalFlags.isSet(CoreGlobalFlags.EDIT_DIST_SINGLE_INDEL_SEEDED_ONLY_FLAG);
  // Enable the heuristic aligners (faster, but some lower quality alignments are produced)
  private static final boolean ENABLE_HEURISTIC_ALIGNING = GlobalFlags.getBooleanValue(CoreGlobalFlags.EDIT_DIST_HEURISTIC_ALIGNERS_FLAG);
  // Specify how many reads to log with -D option
  private static final int EDIT_LOGGING_AMOUNT = GlobalFlags.getIntegerValue(CoreGlobalFlags.EDIT_DIST_LOGGING_AMOUNT_FLAG);

  // Priors for cg CG-gotoh TODO merge with the same default priors that variant calling uses
  private static final String GOTOH_CG_PRIORS = "cg_real_errors";

  /** default table for single indel aligner */
  public static final String DEFAULT_SINGLE_INDEL_TABLE = "alignertable";

  /** The default penalty for a gap open */
  public static final int DEFAULT_GAP_OPEN_PENALTY = 19;
  /** The default penalty for a gap extension */
  public static final int DEFAULT_GAP_EXTEND_PENALTY = 1;
  /** The default penalty for a substitution (mismatch) */
  public static final int DEFAULT_SUBSTITUTION_PENALTY = 9;
  /** The default penalty for unknowns (N nucleotide, aligning off template) */
  public static final int DEFAULT_UNKNOWNS_PENALTY = 5;

  private static final int MIN_LEN_FOR_SEEDED_ALIGNER = 64;

  private EditDistanceFactory() {
  }

  /**
   * Creates a specific Edit Distance implementation depending on the parameters
   *
   * @param ngsParams {@link NgsParams} for current run
   * @param reader1 a sequences reader containing reads
   * @param reader2 another sequences reader containing reads (may be null)
   * @return an edit distance implementation
   */
  public static BidirectionalEditDistance createEditDistance(NgsParams ngsParams, SequencesReader reader1, SequencesReader reader2) {
    if (reader1 == null) {
      throw new IllegalArgumentException("reader1 cannot be null");
    }
    int maxLength = (int) reader1.maxLength();
    int minLength = (int) reader1.minLength();
    final PrereadType prereadType = reader1.getPrereadType();
    if (reader2 != null) {
      if (reader1.getPrereadType() != reader2.getPrereadType()) {
        throw new IllegalArgumentException("Left and right SDF types do not match");
      }
      maxLength = Math.max(maxLength, (int) reader2.maxLength());
      minLength = Math.min(minLength, (int) reader2.minLength());
    }
    if (prereadType == PrereadType.CG) {
      if (reader2 == null) {
        throw new IllegalArgumentException("CG requires paired data");
      }
    }
    return createEditDistance(ngsParams, prereadType, minLength, maxLength);
  }

  /**
   * Creates a specific Edit Distance implementation depending on the parameters
   *
   * @param ngsParams {@link NgsParams} for current run
   * @param prereadType type of the preread
   * @param minReadLength minimum read length
   * @param maxReadLength maximum read length
   * @return an edit distance implementation
   */
  static BidirectionalEditDistance createEditDistance(NgsParams ngsParams, PrereadType prereadType, int minReadLength, int maxReadLength) {

    final NgsOutputParams ngsOutputParams = ngsParams.outputParams();
    final SAMReadGroupRecord samReadGroupRecord = ngsOutputParams == null ? null : ngsOutputParams.readGroup();
    final String platform = samReadGroupRecord == null ? null : samReadGroupRecord.getPlatform();
    final AlignerMode effectiveChain;
    if (ngsParams.alignerMode() == AlignerMode.AUTO) {
      if (MachineType.ILLUMINA_PE.compatiblePlatform(platform) || MachineType.ILLUMINA_SE.compatiblePlatform(platform)) {
        effectiveChain = AlignerMode.TABLE;
      } else {
        effectiveChain = AlignerMode.GENERAL;
      }
      Diagnostic.userLog("Automatically chosen aligner mode: " + effectiveChain.toString().toLowerCase(Locale.getDefault()));
    } else {
      effectiveChain = ngsParams.alignerMode();
    }
    if (USE_SINGLE_INDEL_SEEDED_ONLY) {
      Diagnostic.developerLog("Using SingleIndelSeededEditDistance: maxReadLength=" + maxReadLength);
      return new RcEditDistance(new SoftClipper(new UnidirectionalPrioritisedEditDistance(new SingleIndelSeededEditDistance(ngsParams, maxReadLength)), ngsParams.indelSoftClipDistance(), ngsParams.mismatchSoftClipDistance(), ngsParams.minMatches()));
    } else if (USE_GOTOH_ONLY) {
      Diagnostic.developerLog("Using Gotoh only");
      return new RcEditDistance(new SoftClipper(new GotohEditDistance(ngsParams.gapOpenPenalty(), ngsParams.gapExtendPenalty(), ngsParams.substitutionPenalty(), ngsParams.unknownsPenalty(), false), ngsParams.indelSoftClipDistance(), ngsParams.mismatchSoftClipDistance(), ngsParams.minMatches()));
    } else if (effectiveChain == AlignerMode.TABLE) {
      Diagnostic.developerLog("Using SingleIndelEditDistance (TABLE): maxReadLength=" + maxReadLength);
      return new RcEditDistance(new SoftClipper(new UnidirectionalPrioritisedEditDistance(new SingleIndelEditDistance(ngsParams, maxReadLength)), ngsParams.indelSoftClipDistance(), ngsParams.mismatchSoftClipDistance(), ngsParams.minMatches()));
    }
    // General case

    final boolean useSeededAligner = maxReadLength > MIN_LEN_FOR_SEEDED_ALIGNER;

    if (prereadType == PrereadType.CG) {
      if (minReadLength != maxReadLength) {
        throw new IllegalArgumentException("CG data requires fixed length reads");
      }
      if ((minReadLength == CgUtils.CG_RAW_READ_LENGTH) || (minReadLength == CgUtils.CG2_RAW_READ_LENGTH)) {
        Diagnostic.developerLog("Using Gotoh CG aligner");
        return new RcCgEditDistance(createCgGotohEditDistance(ngsParams.unknownsPenalty(), minReadLength));
      }
      throw new IllegalArgumentException("CG data does not support read length of " + minReadLength);
    } else if (useSeededAligner) {
      Diagnostic.developerLog("Using aligner chain");
      final int lbValue = calculateLowerBoundValue(minReadLength, ngsParams.alignerBandWidthFactor());

      final ArrayList<UnidirectionalEditDistance> fwd = new ArrayList<>();
      final ArrayList<UnidirectionalEditDistance> rev = new ArrayList<>();
      if (ngsParams.outputParams() == null || ngsParams.outputParams().readGroup() == null || !MachineType.IONTORRENT.compatiblePlatform(platform)) {  //if we aren't iontorrent, use noindels edit distance
        Diagnostic.developerLog("NoIndelsEditDistance");
        fwd.add(new NoIndelsEditDistance(ngsParams));
        rev.add(new NoIndelsEditDistance(ngsParams));
      }
      Diagnostic.developerLog("LowerBoundEditDistance lbValue=" + lbValue);
      fwd.add(new LowerBoundEditDistance(lbValue, ngsParams.substitutionPenalty(), ngsParams.unknownsPenalty()));
      rev.add(new LowerBoundEditDistance(lbValue, ngsParams.substitutionPenalty(), ngsParams.unknownsPenalty()));
      Diagnostic.developerLog("LoggingOnlyEditDistance");
      fwd.add(new LoggingOnlyEditDistance(EDIT_LOGGING_AMOUNT));
      rev.add(new LoggingOnlyEditDistance(EDIT_LOGGING_AMOUNT));
      if (ENABLE_HEURISTIC_ALIGNING) {
        Diagnostic.developerLog("HopStepEditDistanceLong");
        fwd.add(new HopStepEditDistanceLong(ngsParams));
        rev.add(new HopStepEditDistanceLong(ngsParams));
      }
      Diagnostic.developerLog("SeededAligner");
      fwd.add(new SeededAligner(ngsParams, !ENABLE_HEURISTIC_ALIGNING));
      rev.add(new SeededAligner(ngsParams, !ENABLE_HEURISTIC_ALIGNING));
      Diagnostic.developerLog("GotohEditDistance");
      fwd.add(new GotohEditDistance(ngsParams.gapOpenPenalty(), ngsParams.gapExtendPenalty(), ngsParams.substitutionPenalty(), ngsParams.unknownsPenalty(), false));
      rev.add(new GotohEditDistance(ngsParams.gapOpenPenalty(), ngsParams.gapExtendPenalty(), ngsParams.substitutionPenalty(), ngsParams.unknownsPenalty(), false));

      return new RcEditDistance(
        new SoftClipper(new UnidirectionalPrioritisedEditDistance(fwd.toArray(new UnidirectionalEditDistance[0])),
          ngsParams.indelSoftClipDistance(), ngsParams.mismatchSoftClipDistance(), ngsParams.minMatches()),
        new SoftClipper(new UnidirectionalPrioritisedEditDistance(rev.toArray(new UnidirectionalEditDistance[0])),
          ngsParams.indelSoftClipDistance(), ngsParams.mismatchSoftClipDistance(), ngsParams.minMatches()));
    }

    return new RcEditDistance(new SoftClipper(new UnidirectionalPrioritisedEditDistance(
      new NoIndelsEditDistance(ngsParams),
      new GotohEditDistance(ngsParams.gapOpenPenalty(), ngsParams.gapExtendPenalty(), ngsParams.substitutionPenalty(), ngsParams.unknownsPenalty(), false)),
      ngsParams.indelSoftClipDistance(), ngsParams.mismatchSoftClipDistance(), ngsParams.minMatches()));
  }

  private static CgGotohEditDistance createCgGotohEditDistance(int unknownsPenalty, int readLength) {
    try {
      final boolean v2 = readLength == CgUtils.CG2_RAW_READ_LENGTH;
      return new CgGotohEditDistance(7, new RealignParamsImplementation(new MachineErrorParamsBuilder().errors(GOTOH_CG_PRIORS).create()), unknownsPenalty, v2);
    } catch (final InvalidParamsException | IOException e) {
      throw new RuntimeException(e);
    }
  }

  /**
   * Calculate how wide the lower-bound window should be.
   *
   * @param minReadLength minimum read length
   * @param maxShiftFactor determines max shift for the given read length
   * @return width of lower bound window (word size).
   */
  protected static int calculateLowerBoundValue(int minReadLength, MaxShiftFactor maxShiftFactor) {
    //    final String lowerBoundValue = System.getProperty(LB_VALUE);

    // auto selection of the lower bound value
    int lbValue = 3;
    final StringBuilder sb = new StringBuilder();
    sb.append("EditDistanceFactory: auto lower bound determination");

    final int minscore = 1 + 2 * maxShiftFactor.calculateMaxShift(minReadLength);
    //    final int minscore = 1 + 2 * Math.min(MaxShiftUtils.calculateDefaultMaxShift(minReadLength), maxScore); // 100 long 10% = 10

    for (int width = 3; width < 8; ++width) {
      final long hashtotal = (long) Math.pow(4, width); // e.g. for width=3 it's hashtotal == 64
      final long numhashes = minReadLength / width; // e.g. for 100 long and width = 3, then that's 33 hashes

      sb.append("EditDistanceFactory: lowerbound width=").append(width).append(", hashtotal=").append(hashtotal)
      .append(" numhashes=").append(numhashes).append(" minscore=").append(minscore).append(StringUtils.LS);
      // the total number of potential hashes should be much greater than the minscore
      // and the number of hashes that are calculated in the read need to be greater than minscore too.
      if ((hashtotal > 6L * minscore) && (numhashes > minscore)) {
        lbValue = width;
        break;
      }
    }
    sb.append("EditDistanceFactory: selecting lbValue=").append(lbValue);

    // override the lower bound using a -D option
    //    if (lowerBoundValue != null) {
    //      lbValue = Integer.parseInt(lowerBoundValue);
    //      sb.append("EditDistanceFactory: overridden by -D option to ").append(lbValue);
    //    }
    Diagnostic.developerLog(sb.toString());
    return lbValue;
  }

  /**
   * Factory method for creating a protein edit distance object.
   * @param matrix one of the protein scoring matrices.
   * @return a GotohProteinEditDistance.
   */
  public static BidirectionalEditDistance createProteinEditDistance(final ProteinScoringMatrix matrix) {
    return new UnidirectionalAdaptor(new GotohProteinEditDistance(matrix));
  }
}
