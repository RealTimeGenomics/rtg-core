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
