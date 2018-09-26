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

  enum Classifications {
    /** Treat all instances in this dataset as positive examples */
    ALL_POSITIVE,
    /** Treat all instances in this dataset as negative examples */
    ALL_NEGATIVE,
    /** Determine class label for each instance based on INFO annotations */
    ANNOTATED
  }

  private final File mVcfFile;
  private final int mSampleNum;
  private final Classifications mClassType;
  private final double mInstanceWeight;
  private final boolean mReweight;

  /**
   * Constructor
   * @param vcfFile the file containing training variants
   * @param sampleNum the index of the sample from which to obtain format-level attributes
   * @param classType specifies how to determine the classification for instances in this dataset
   * @param reweight true if positive and negative instances should be equalized in weight
   * @param instanceWeight the weight to assign to instances from this dataset.
   */
  public VcfDataset(File vcfFile, int sampleNum, Classifications classType, boolean reweight, double instanceWeight) {
    mVcfFile = vcfFile;
    mSampleNum = sampleNum;
    mClassType = classType;
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
    return mClassType == Classifications.ALL_POSITIVE;
  }

  /**
   * @return the method that should be used to determine the class label of each instance in the dataset.
   */
  public Classifications classifications() {
    return mClassType;
  }

  public boolean isReweight() {
    return mReweight;
  }

  public double getInstanceWeight() {
    return mInstanceWeight;
  }

  @Override
  public String toString() {
    return "VcfDataset "
      + mClassType
      + " reweight=" + mReweight
      + " sample=" + mSampleNum
      + " weight=" + mInstanceWeight
      + " " + mVcfFile;
  }
}
