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
