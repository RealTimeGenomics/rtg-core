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

import com.rtg.variant.bayes.CodeDiploid;
import com.rtg.variant.bayes.CodeHaploid;
import com.rtg.variant.bayes.Evidence;
import com.rtg.variant.bayes.EvidenceInterface;

import junit.framework.TestCase;

/**
 */
public class EvidenceQTest extends TestCase {

  protected static void check(final EvidenceInterface p, final double expected, final int hy) {
    assertEquals("expected=" + expected + " hy=" + hy, expected, p.probability(hy), 0.00001);
  }

  public void testDiploidEqualsHaploid() {
    // Evidence Q relies on the property that the homozygous portion of the diploid code is indexed identically to the haploid code
    final int max = 7;
    final CodeDiploid diploid = new CodeDiploid(max);
    final CodeHaploid haploid = new CodeHaploid(max);
    for (int i = 0; i < max; ++i) {
      assertTrue(diploid.homozygous(i));
      assertEquals(haploid.a(i), diploid.a(i));
    }

  }

  public void test() {
    final Evidence dq = new EvidenceQ(DescriptionSnp.SINGLETON, 1, 0, 0, 0.3, 0.1, true, false, true, false, false);
    dq.integrity();

    check(dq, 0.1 / 3.0, 0);
    check(dq, 1.0 - 0.1, 1);
    check(dq, 0.1 / 3.0, 2);
    check(dq, 0.1 / 3.0, 3);
    for (int i = 0; i < 4; ++i) {
      assertEquals(0.25, dq.pe());
    }
    assertEquals(0.1, dq.error());
    assertEquals(1, dq.read());
    assertFalse(dq.isUnmapped());
  }

  public void testUnmapped() {
    final Evidence dq = new EvidenceQ(DescriptionSnp.SINGLETON, 1, 0, 0, 0.3, 0.1, true, false, true, false, true);
    assertTrue(dq.isUnmapped());
  }
}
