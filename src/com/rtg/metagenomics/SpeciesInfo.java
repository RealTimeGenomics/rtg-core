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

package com.rtg.metagenomics;

/**
 * Information about a species.
 */
public class SpeciesInfo {
  final double mCoverageDepth;
  final long mCoverageBreadth;
  final double mMappedReads;

  /**
   * @return Returns the coverage Depth.
   */
  public double getCoverageDepth() {
    return mCoverageDepth;
  }

  /**
   * @return Returns the coverage Breadth.
   */
  public long getCoverageBreadth() {
    return mCoverageBreadth;
  }

  /**
   * @return Returns the mapped Reads.
   */
  public double getMappedReads() {
    return mMappedReads;
  }

  /**
   * @param coverageDepth coverage depth.
   * @param coveragedBreadth  coverage depth.
   * @param mappedReads number of mapped reads.
   */
  public SpeciesInfo(double coverageDepth, long coveragedBreadth, double mappedReads) {
    mCoverageDepth = coverageDepth;
    mCoverageBreadth = coveragedBreadth;
    mMappedReads = mappedReads;
  }
}
