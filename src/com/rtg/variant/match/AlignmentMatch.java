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
package com.rtg.variant.match;

import com.rtg.ngs.Arm;
import com.rtg.reader.FastaUtils;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.Integrity;
import com.rtg.variant.MachineErrorChooserInterface;
import com.rtg.variant.NoQualityException;
import com.rtg.variant.PhredScaler;
import com.rtg.variant.VariantAlignmentRecord;
import com.rtg.variant.dna.DNARange;
import com.rtg.variant.dna.DNARangeNAT;
import com.rtg.variant.util.VariantUtils;

/**
 * Represent a match between a read (alignment) and reference.
 */
public class AlignmentMatch extends Match implements Integrity {

  private static final DNARange DNA_RANGE = DNARange.RANGE;

  private final double mMapError;
  private final int[] mReadNt;
  private final double[] mBaseError;
  private final double mQualityDefault;
  private final String mReadString;
  private final boolean mFixedLeft;
  private final boolean mFixedRight;
  private int mBasesLeftOfMatch;
  private int mBasesRightOfMatch;
  private final VariantAlignmentRecord mAlignmentRecord;


  /**
   * For testing only, as it uses ascii qualities and does not use a <code>MachineErrorChooser</code>, so does no quality curve correction.
   * A new single nucleotide difference match.
   *
   * @param var associated alignment record (if any).
   * @param readNt read nucleotides
   * @param qScore Phred quality of read nucleotides
   * @param qDefault default Phred quality of read nucleotides Used if <code>qScore</code> null)
   * @param start zero-based offset into read nucleotides
   * @param length length of read region
   * @param readScore Phred quality of read mapping
   */
  public AlignmentMatch(final VariantAlignmentRecord var, final String readNt, final String qScore, final int qDefault, final int start, final int length, final int readScore) {
    this(var, null, readNt, FastaUtils.asciiToRawQuality(qScore), qDefault, start, length, readScore, true, true);
  }

  /**
   * A new single nucleotide difference match
   * @param var associated alignment record (if any).
   * @param chooser Machine error chooser (if null, no quality correction is done)
   * @param readNt read nucleotides
   * @param qScore binary Phred quality of read nucleotides
   * @param qDefault default Phred quality of read nucleotides (used if <code>qScore</code> null)
   * @param start zero-based offset into read nucleotides
   * @param length length of read region
   * @param readScore Phred quality of read mapping
   * @param fixedLeft true if read is fixed (anchored) at left end.
   * @param fixedRight true if read is fixed (anchored) at right end.
   */
  public AlignmentMatch(final VariantAlignmentRecord var, final MachineErrorChooserInterface chooser, final String readNt, final byte[] qScore, final int qDefault, final int start, final int length, final int readScore, final boolean fixedLeft, final boolean fixedRight) {
    super();
    mAlignmentRecord = var;
    final PhredScaler me = chooser == null ? null : chooser.machineErrors(var.getReadGroup(), var.isReadPaired());
    mMapError = VariantUtils.phredToProb(readScore);
    if (qScore == null) {
      mBaseError = null;
    } else {
      mBaseError = new double[length];
      for (int k = 0; k < mBaseError.length; k++) {
        // apply machine error calibration curve for appropriate read group if possible
        final int phred = qScore[k + start];
        assert 0 <= phred;
        mBaseError[k] = VariantUtils.phredToProb(phred);
      }
    }
    // apply correct machine error calibration curve to qDefault, if possible
    final Arm arm  = Arm.LEFT;
    mQualityDefault = VariantUtils.phredToProb(me == null ? qDefault : me.getScaledPhred((byte) qDefault, 0, arm));
    mFixedLeft = fixedLeft;
    mFixedRight = fixedRight;
    mReadNt = new int[length];
    final StringBuilder sb = new StringBuilder();
    for (int k = 0; k < length; k++) {
      final int v = DNA_RANGE.valueOf(readNt.charAt(k + start));
      mReadNt[k] = v;
      sb.append(DNA_RANGE.toString(v));
    }
    mReadString = sb.toString();
  }

  /**
   * Get the associated SAM record.
   * @return the associated SAM record (may be null for legacy reasons).
   */
  public VariantAlignmentRecord alignmentRecord() {
    return mAlignmentRecord;
  }


  @Override
  public double baseError(final int index) throws NoQualityException, IndexOutOfBoundsException {
    if (mBaseError == null) {
      return mQualityDefault;
    }
    return mBaseError[index];
  }

  @Override
  public int read(final int index) {
    final int res = mReadNt[index] - 1;
    assert DNARangeNAT.DNA.valid(res);
    return res;
  }

  public int getBasesLeftOfMatch() {
    return mBasesLeftOfMatch;
  }

  public void setBasesLeftOfMatch(int basesLeftOfMatch) {
    mBasesLeftOfMatch = basesLeftOfMatch;
  }

  public int getBasesRightOfMatch() {
    return mBasesRightOfMatch;
  }

  public void setBasesRightOfMatch(int basesRightOfMatch) {
    mBasesRightOfMatch = basesRightOfMatch;
  }

  @Override
  public double mapError() {
    return mMapError;
  }

  @Override
  public int length() {
    return mReadNt.length;
  }

  @Override
  public String readString() {
    return mReadString;
  }

  @Override
  public boolean isFixedLeft() {
    return mFixedLeft;
  }

  @Override
  public boolean isFixedRight() {
    return mFixedRight;
  }

  @Override
  public boolean integrity() {
    Exam.assertEquals(mReadString, mReadString.length(), mReadNt.length);
    if (mBaseError != null) {
      Exam.assertEquals(mReadNt.length, mBaseError.length);
      Exam.assertProbabilities(mBaseError);
    }
    Exam.assertProbability(mMapError);
    Exam.assertProbability(mQualityDefault);
    return true;
  }

  @Override
  public boolean globalIntegrity() {
    return true;
  }
}
