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
package com.rtg.variant.bayes.snp;

import com.rtg.variant.bayes.EvidenceInterface;

import junit.framework.TestCase;

/**
 */
public class IndelDetectorTest extends TestCase {

  //single identity insertion
  public void testa() {
    final IndelDetector dis = new IndelDetector();

    final EvidenceInterface prob = new EvidenceIndel(0.42, 0, 0);
    dis.increment(prob);

    assertEquals(0, dis.nonTrivialDeletionCount());
    assertEquals(1, dis.nonTrivialInsertCount());
  }

  public void test() throws Exception {
    check(0);
    check(1);
    check(2);
    check(10);
    check(20);
    check(40);
  }

  private void check(final int nonTrivial) {
    final IndelDetector dis = new IndelDetector();
    for (int i = 0; i < nonTrivial; ++i) {
      dis.increment(new EvidenceIndel(0.0, 0, 0));
    }

    assertEquals(nonTrivial, dis.nonTrivialInsertCount());
  }

  public void testIndelCount() {
    final IndelDetector dis = new IndelDetector();
    final double fraction = 0.05;

    for (int i = 0; i < 10; ++i) {
      dis.increment(null);
      dis.increment(new EvidenceIndel(0.0, 0, 0));
    }
    assertEquals(2, dis.minIndelCount(fraction)); //10

    for (int i = 0; i < 89; ++i) {
      dis.increment(null);
      dis.increment(new EvidenceIndel(0.0, 0, 0));
    }
    assertEquals(4, dis.minIndelCount(fraction));  //99

    dis.increment(null);
    assertEquals(2, dis.minIndelCount(0));  //100
    assertEquals(5, dis.minIndelCount(fraction));

    for (int i = 0; i < 72; ++i) {
      dis.increment(null);
      dis.increment(new EvidenceIndel(0.0, 0, 0));
    }
    assertEquals(8, dis.minIndelCount(fraction)); //172
  }
}
