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
package com.rtg.variant.util;

import java.io.IOException;

import junit.framework.TestCase;

/**
 */
public class TrioConcordanceTest extends TestCase {

  public void testConcordance() throws IOException {
    TrioConcordance tc = new TrioConcordance("child", "father", "mother");

    for (int i = 0; i < 25; i++) {
      tc.add(new Genotype("0/0"), new Genotype("0/0"), new Genotype("0/0"));
      tc.add(new Genotype("1/0"), new Genotype("0/1"), new Genotype("1/0"));
      tc.add(new Genotype("1/0"), new Genotype("0/1"), new Genotype("0/0"));
      tc.add(new Genotype("1/1"), new Genotype("0/1"), new Genotype("1/1"));
    }

    assertEquals(TrioConcordance.Status.OK, tc.getStatus(100, 99.0));

    tc.add(new Genotype("1/1"), new Genotype("0/0"), new Genotype("1/1"));
    assertEquals(TrioConcordance.Status.MOTHER, tc.getStatus(100, 99.5));
    assertEquals(TrioConcordance.Status.OK, tc.getStatus(200, 99.5));

    tc.add(new Genotype("1/1"), new Genotype("0/0"), new Genotype("0/0"));
    assertEquals(TrioConcordance.Status.BOTH, tc.getStatus(100, 99.5));
    assertEquals(TrioConcordance.Status.OK, tc.getStatus(200, 99.5));

    tc.add(new Genotype("1/1"), new Genotype("0/."), new Genotype("0/0"));
    tc.add(new Genotype("1/1"), new Genotype("./."), new Genotype("0/0"));
    assertEquals(TrioConcordance.Status.OK, tc.getStatus(105, 99.5));
    tc.add(new Genotype("1/1"), new Genotype("0/."), new Genotype("0/0"));
    assertEquals(TrioConcordance.Status.FATHER, tc.getStatus(105, 99.5));
  }

}
