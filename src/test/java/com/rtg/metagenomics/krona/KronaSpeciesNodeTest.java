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

import junit.framework.TestCase;

/**
 */
public class KronaSpeciesNodeTest extends TestCase {
  public void testIt() {
    final BoundedDouble ab = new BoundedDouble(0.5, 0.2, 0.6);
    final BoundedDouble dna = new BoundedDouble(0.2, 0.1, 0.4);
    final KronaSpeciesNode node = new KronaSpeciesNode(ab, dna, 2.0, 40.2, 2.3, 5.0, (long) 7);

    assertEquals(2.0, node.mConfidence);
    assertEquals(40.2, node.mMappedReads);

    assertEquals(ab, node.mAbundance);
    assertEquals(dna, node.mDnaFraction);
    assertEquals(2.3, node.mCoverageDepth);
    assertEquals(5.0, node.mCoverageBreadth);
    assertEquals(Long.valueOf(7), node.mGenomeLength);
  }
}
