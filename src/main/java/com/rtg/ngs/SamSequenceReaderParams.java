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
package com.rtg.ngs;

/**
 * Parameters for loading sequence data from sam files
 */
public class SamSequenceReaderParams {
  private final boolean mFlattenPairs;
  private final boolean mUnorderedLoad;

  /**
   * @param flattenPairs whether when running in single end mode to accept paired end data and load both sides into single SDF
   * @param unorderedLoad support loading SAM files that are not read name sorted at the cost of memory
   */
  public SamSequenceReaderParams(boolean flattenPairs, boolean unorderedLoad) {
    mFlattenPairs = flattenPairs;
    mUnorderedLoad = unorderedLoad;
  }

  /**
   * @return whether when running in single end mode to accept paired end data and load both sides into single SDF
   */
  public boolean flattenPairs() {
    return mFlattenPairs;
  }

  /**
   * @return support loading SAM files that are not read name sorted at the cost of memory
   */
  public boolean unorderedLoad() {
    return mUnorderedLoad;
  }

}
