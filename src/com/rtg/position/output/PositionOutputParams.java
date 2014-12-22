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
package com.rtg.position.output;

import java.io.File;

import com.rtg.launcher.OutputDirParams;
import com.rtg.util.ObjectParams;
import com.rtg.util.Utils;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.Integrity;

/**
 * Specifies output parameters for Position.
 */
public class PositionOutputParams extends ObjectParams implements OutputDirParams, Integrity {

  private final File mOutDir;

  private final OutputFormatType mFormatType;

  private final PositionDistributionParams mDistrParams;

  private final Double mScoreThreshold;

  private final boolean mZip;

  private final int mTopN;

  /**
   * @param out where output is to be written.
   * @param type specifies format of output.
   * @param distr parameters for specifying the probability distribution used in gapped extension calculations.
   * @param scoreThreshold any SNPs with score less than this are ignored.
   * @param zip true iff output to be zipped.
   * @param topN maximum results to retain per query sequence.
   */
  public PositionOutputParams(final File out, final OutputFormatType type, final PositionDistributionParams distr, final Double scoreThreshold, final boolean zip, final int topN) {
    mOutDir = out;
    mFormatType = type;
    mDistrParams = distr;
    mScoreThreshold = scoreThreshold;
    mZip = zip;
    mTopN = topN;
    mObjects = new Object[] {mOutDir, mFormatType, mDistrParams, mScoreThreshold, mZip, mTopN};
  }

  /**
   * Get where the output is to go.
   * @return where the output is to go.
   */
  @Override
  public File directory() {
    return mOutDir;
  }

  /**
   * Get the output format.
   * @return the output format.
   */
  public OutputFormatType format() {
    return mFormatType;
  }

  /**
   * Get the maximum allowed gap.
   * @return the maximum allowed gap.
   */
  public Integer maxGap() {
    if (mDistrParams == null) {
      return null;
    }
    return mDistrParams.maxGap();
  }

  /**
   * Get the distribution parameters.
   * @return the distribution parameters.
   */
  public PositionDistributionParams distribution() {
    return mDistrParams;
  }

  /**
   * Get the threshold for scores (used in similarity and region).
   *
   * @return the threshold for scores.
   */
  public Double scoreThreshold() {
    return mScoreThreshold;
  }

  /**
   * Test if output should be zipped.
   *
   * @return true iff output should be zipped.
   */
  public boolean zip() {
    return mZip;
  }

  /**
   * Return the maximum number of results per query to retain.
   *
   * @return maximum results
   */
  public int topN() {
    return mTopN;
  }

  /**
   * Check if a and b are equal when either can be null.
   * @return true iff a and b are equal including being both null.
   */
  static boolean equalNull(final Object a, final Object b) {
    if (a == b) {
      return true;
    }
    if (a == null || b == null) {
      return false;
    }
    return a.equals(b);
  }

  @Override
  public String toString() {
    return "outputDir=" + mOutDir + " format=" + mFormatType + " zip=" + mZip
      + " score threshold=" + Utils.realFormat(mScoreThreshold) + " topN=" + mTopN
      + " distribution={" + Utils.nullFormat(mDistrParams) + "}";
  }

  @Override
  public boolean integrity() {
    Exam.assertTrue(mOutDir != null);
    Exam.assertNotNull(mFormatType);
    return true;
  }

  @Override
  public boolean globalIntegrity() {
    integrity();
    return true;
  }
}
