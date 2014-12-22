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
package com.rtg.metagenomics.krona;

import com.rtg.util.BoundedDouble;

/**
 * Stores information necessary to construct a node for the species <code>Krona</code> report
 */
public class KronaSpeciesNode {
  final BoundedDouble mAbundance;
  final BoundedDouble mDnaFraction;

  final Double mConfidence;
  final Double mCoverageDepth;
  final Double mCoverageBreadth;
  final Double mMappedReads;
  final Long mGenomeLength;

  /**
   * Create a node which stores data sufficient for a <code>Krona</code> element
   * @param abundance the abundance values
   * @param dnaFraction the DNA fraction values
   * @param confidence the confidence
   * @param mappedReads the number of mapped reads
   * @param coverageDepth the coverage depth
   * @param coverageBreadth the coverage breadth
   * @param genomeLength the length of the species genome
   */
  public KronaSpeciesNode(BoundedDouble abundance, BoundedDouble dnaFraction, Double confidence, Double mappedReads, Double coverageDepth, Double coverageBreadth, Long genomeLength) {
    mAbundance = abundance;
    mDnaFraction = dnaFraction;
    mConfidence = confidence;
    mMappedReads = mappedReads;
    mCoverageBreadth = coverageBreadth;
    mCoverageDepth = coverageDepth;
    mGenomeLength = genomeLength;
  }
}
