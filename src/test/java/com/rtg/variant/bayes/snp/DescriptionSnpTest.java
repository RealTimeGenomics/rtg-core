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

import com.rtg.variant.bayes.Description;

import junit.framework.TestCase;

/**
 */
public class DescriptionSnpTest extends TestCase {

  private void check(final String exp, final int index, Description description) {
    assertTrue(description.valid(index));
    assertEquals(exp, description.name(index));
    final StringBuilder sb = new StringBuilder();
    description.writeName(sb, index);
    assertEquals(exp, sb.toString());
  }

  public void test() {
    assertEquals(4, DescriptionSnp.SINGLETON.size());
    assertEquals(1, DescriptionSnp.SINGLETON.minLength());
    assertEquals(1, DescriptionSnp.SINGLETON.maxLength());
    check("A", 0, DescriptionSnp.SINGLETON);
    check("C", 1, DescriptionSnp.SINGLETON);
    check("G", 2, DescriptionSnp.SINGLETON);
    check("T", 3, DescriptionSnp.SINGLETON);

    assertFalse(DescriptionSnp.SINGLETON.valid(-1));
    assertFalse(DescriptionSnp.SINGLETON.valid(4));
  }
}
