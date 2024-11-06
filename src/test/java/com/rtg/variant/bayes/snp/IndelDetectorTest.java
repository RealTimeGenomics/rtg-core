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

  public void test() {
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
