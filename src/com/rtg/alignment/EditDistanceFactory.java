/*
 * Copyright (c) 2014. Real Time Genomics Limited.
 *
 * Use of this source code is bound by the Real Time Genomics Limited Software Licence Agreement
 * for Academic Non-commercial Research Purposes only.
 *
 * If you did not receive a license accompanying this file, a copy must first be obtained by email
 * from support@realtimegenomics.com.  On downloading, using and/or continuing to use this source
 * code you accept the terms of that license agreement and any amendments to those terms that may
 * be made from time to time by Real Time Genomics Limited.
 */
package com.rtg.alignment;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Locale;

import com.rtg.launcher.GlobalFlags;
import com.rtg.mode.ProteinScoringMatrix;
import com.rtg.ngs.NgsOutputParams;
import com.rtg.ngs.NgsParams;
import com.rtg.protein.GotohProteinEditDistance;
import com.rtg.reader.PrereadType;
import com.rtg.reader.SequencesReader;
import com.rtg.util.InvalidParamsException;
import com.rtg.util.MaxShiftFactor;
import com.rtg.util.StringUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.machine.MachineType;
import com.rtg.variant.MachineErrorParamsBuilder;
import com.rtg.variant.realign.RealignParamsImplementation;

import net.sf.samtools.SAMReadGroupRecord;

/**
 * Factory which creates a specific Edit Distance implementation depending on
 * the parameters passed.
 *
 */
public final class EditDistanceFactory {

  //  private static final String LB_VALUE = "lower-bound-value";

  // Only use the gotoh aligner (disable all others, has priority over USE_NOINDELS_ONLY)
  private static final boolean USE_GOTOH_ONLY = GlobalFlags.isSet(GlobalFlags.EDIT_DIST_GOTOH_ONLY_FLAG);
  // Only use the single-indel-seeded aligner (disable all others, has priority over USE_NOINDELS_ONLY)
  private static final boolean USE_SINGLE_INDEL_SEEDED_ONLY = GlobalFlags.isSet(GlobalFlags.EDIT_DIST_SINGLE_INDEL_SEEDED_ONLY_FLAG);
  // Enable the heuristic aligners (faster, but some lower quality alignments are produced)
  private static final boolean ENABLE_HEURISTIC_ALIGNING = GlobalFlags.getBooleanValue(GlobalFlags.EDIT_DIST_HEURISTIC_ALIGNERS_FLAG);
  // Specify how many reads to log with -D option
  private static final int EDIT_LOGGING_AMOUNT = GlobalFlags.getIntegerValue(GlobalFlags.EDIT_DIST_LOGGING_AMOUNT_FLAG);

  // Priors for cg CG-gotoh
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
  public static EditDistance createEditDistance(NgsParams ngsParams, SequencesReader reader1, SequencesReader reader2) {
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
  static EditDistance createEditDistance(NgsParams ngsParams, PrereadType prereadType, int minReadLength, int maxReadLength) {

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
      return new SoftClipperOmni(new RcEditDistance(new UnidirectionalPrioritisedEditDistance(ngsParams, new SingleIndelSeededEditDistance(ngsParams, maxReadLength))), ngsParams.softClipDistance());
    } else if (USE_GOTOH_ONLY) {
      Diagnostic.developerLog("Using Gotoh only");
      return new SoftClipperOmni(new RcEditDistance(new GotohEditDistance(ngsParams)), ngsParams.softClipDistance());
    } else if (effectiveChain == AlignerMode.TABLE) {
      Diagnostic.developerLog("Using SingleIndelEditDistance (TABLE): maxReadLength=" + maxReadLength);
      return new SoftClipperOmni(new RcEditDistance(new UnidirectionalPrioritisedEditDistance(ngsParams, new SingleIndelEditDistance(ngsParams, maxReadLength))), ngsParams.softClipDistance());
    }
    // General case

    final boolean useSeededAligner = maxReadLength > MIN_LEN_FOR_SEEDED_ALIGNER;

    if (prereadType == PrereadType.CG) {
      if ((maxReadLength == CgGotohEditDistance.CG_RAW_READ_LENGTH) && (minReadLength == maxReadLength)) {
        Diagnostic.developerLog("Using Gotoh CG aligner");
        return new RcEditDistance(createCgGotohEditDistance(ngsParams.unknownsPenalty()));
      }
      throw new IllegalArgumentException("CG data only supports " + CgGotohEditDistance.CG_RAW_READ_LENGTH + " length reads");
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
      fwd.add(new GotohEditDistance(ngsParams));
      rev.add(new GotohEditDistance(ngsParams));

      return new SoftClipperOmni(new RcEditDistance(
          new UnidirectionalPrioritisedEditDistance(ngsParams, fwd.toArray(new UnidirectionalEditDistance[fwd.size()])),
          new UnidirectionalPrioritisedEditDistance(ngsParams, rev.toArray(new UnidirectionalEditDistance[rev.size()]))), ngsParams.softClipDistance());
    }

    return new SoftClipperOmni(new RcEditDistance(new UnidirectionalPrioritisedEditDistance(ngsParams,
        new NoIndelsEditDistance(ngsParams),
        new GotohEditDistance(ngsParams))), ngsParams.softClipDistance());
  }

  private static CgGotohEditDistance createCgGotohEditDistance(int unknownsPenalty) {
    try {
      return new CgGotohEditDistance(7, new RealignParamsImplementation(new MachineErrorParamsBuilder().errors(GOTOH_CG_PRIORS).create()), unknownsPenalty);
    } catch (final InvalidParamsException | IOException e) {
      throw new RuntimeException(e);
    }
  }

  /**
   * Calculate how wide the lower-bound window should be.
   *
   * @param minReadLength minimum read length
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

    for (int width = 3; width < 8; width++) {
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
  public static EditDistance createProteinEditDistance(final ProteinScoringMatrix matrix) {
    return new GotohProteinEditDistance(matrix);
  }
}
