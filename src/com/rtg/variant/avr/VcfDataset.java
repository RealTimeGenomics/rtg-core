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
package com.rtg.variant.avr;

import java.io.File;

/**
 * Encapsulates a single-VCF training dataset
 */
public class VcfDataset {

  private final File mVcfFile;
  private final int mSampleNum;
  private final boolean mIsPositive;
  private final double mInstanceWeight;
  private final boolean mReweight;

  /**
   * Constructor
   * @param vcfFile the file containing training variants
   * @param sampleNum the index of the sample from which to obtain format-level attributes
   * @param isPositive true if the dataset contains positive examples
   * @param reweight true if positive and negative instances should be equalized in weight
   * @param instanceWeight the weight to assign to instances from this dataset.
   */
  public VcfDataset(File vcfFile, int sampleNum, boolean isPositive, boolean reweight, double instanceWeight) {
    mVcfFile = vcfFile;
    mSampleNum = sampleNum;
    mIsPositive = isPositive;
    mReweight = reweight;
    mInstanceWeight = instanceWeight;
  }

  public File getVcfFile() {
    return mVcfFile;
  }

  public int getSampleNum() {
    return mSampleNum;
  }

  public boolean isPositive() {
    return mIsPositive;
  }

  public boolean isReweight() {
    return mReweight;
  }

  public double getInstanceWeight() {
    return mInstanceWeight;
  }
}
