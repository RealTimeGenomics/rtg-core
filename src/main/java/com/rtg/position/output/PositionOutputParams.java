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
