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
package com.rtg.variant.avr;

import java.io.File;

import com.rtg.util.intervals.ReferenceRanges;

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
  private ReferenceRanges<String> mRanges;

  /**
   * Constructor
   * @param vcfFile the file containing training variants
   * @param sampleNum the index of the sample from which to obtain format-level attributes
   * @param classType specifies how to determine the classification for instances in this dataset
   * @param reweight true if positive and negative instances should be equalized in weight
   * @param instanceWeight the weight to assign to instances from this dataset.
   */
  public VcfDataset(File vcfFile, int sampleNum, Classifications classType, boolean reweight, double instanceWeight) {
    this(vcfFile, sampleNum, classType, reweight, instanceWeight, null);
  }

  /**
   * Constructor
   * @param vcfFile the file containing training variants
   * @param sampleNum the index of the sample from which to obtain format-level attributes
   * @param classType specifies how to determine the classification for instances in this dataset
   * @param reweight true if positive and negative instances should be equalized in weight
   * @param instanceWeight the weight to assign to instances from this dataset.
   * @param ranges regions to load training variants from, or null for no restriction
   */
  public VcfDataset(File vcfFile, int sampleNum, Classifications classType, boolean reweight, double instanceWeight, ReferenceRanges<String> ranges) {
    mVcfFile = vcfFile;
    mRanges = ranges;
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

  /**
   * @return the training data regions, if any
   */
  public ReferenceRanges<String> ranges() {
    return mRanges;
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
