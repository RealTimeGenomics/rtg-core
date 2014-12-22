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
    for (int i = 0; i < max; i++) {
      assertTrue(diploid.homozygous(i));
      assertEquals(haploid.a(i), diploid.a(i));
    }

  }

  public void test() {
    final Evidence dq = new EvidenceQ(DescriptionSnp.SINGLETON, 1, 0, 0, 0.3, 0.1, true, false, false, false);
    dq.integrity();

    check(dq, 0.1 / 3.0, 0);
    check(dq, 1.0 - 0.1, 1);
    check(dq, 0.1 / 3.0, 2);
    check(dq, 0.1 / 3.0, 3);
    for (int i = 0; i < 4; i++) {
      assertEquals(0.25, dq.pe());
    }
    assertEquals(0.1, dq.error());
    assertEquals(1, dq.read());
    assertFalse(dq.isUnmapped());
  }

  public void testUnmapped() {
    final Evidence dq = new EvidenceQ(DescriptionSnp.SINGLETON, 1, 0, 0, 0.3, 0.1, true, false, false, true);
    assertTrue(dq.isUnmapped());
  }
}
