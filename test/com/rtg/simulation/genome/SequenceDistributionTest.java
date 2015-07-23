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

package com.rtg.simulation.genome;

import java.io.IOException;

import com.rtg.reader.ReaderTestUtils;
import com.rtg.reader.SequencesReader;

import junit.framework.TestCase;

/**
 */
public class SequenceDistributionTest extends TestCase {
  static final double E = 0.1E-9;
  static final double[] PROB = {0.1, 0.3, 0.1, 0.2, 0.2, 0.1};
  public void testDistribution() {
    final SequenceDistribution dist = new SequenceDistribution(PROB);
    check(dist);
  }
  private void check(SequenceDistribution dist) {
    assertEquals(0, dist.selectSequence(0.1 - E));
    assertEquals(1, dist.selectSequence(0.4 - E));
    assertEquals(2, dist.selectSequence(0.5 - E));
    assertEquals(3, dist.selectSequence(0.7 - E));
    assertEquals(4, dist.selectSequence(0.9 - E));
    assertEquals(5, dist.selectSequence(1.0 - E));

    assertEquals(1, dist.selectSequence(0.1 + E));
    assertEquals(2, dist.selectSequence(0.4 + E));
    assertEquals(3, dist.selectSequence(0.5 + E));
    assertEquals(4, dist.selectSequence(0.7 + E));
    assertEquals(5, dist.selectSequence(0.9 + E));

  }

  public void testDefaults() throws IOException {
    final SequencesReader readerDnaMemory = ReaderTestUtils.getReaderDnaMemory(ReaderTestUtils.fasta("A", "AAA", "A", "AA", "AA", "A"));
    check(SequenceDistribution.createDistribution(readerDnaMemory, null));
  }

  public void testOverrideDefaults() throws IOException {
    final SequencesReader readerDnaMemory = ReaderTestUtils.getReaderDnaMemory(ReaderTestUtils.fasta("AA", "AA", "AA", "AA", "AA", "AA"));
    check(SequenceDistribution.createDistribution(readerDnaMemory, PROB));
  }
}
