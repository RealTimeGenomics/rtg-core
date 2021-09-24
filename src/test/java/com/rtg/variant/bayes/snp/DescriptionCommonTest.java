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

import junit.framework.TestCase;

/**
 */
public class DescriptionCommonTest extends TestCase {

  public void test() {
    final DescriptionCommon de = new DescriptionCommon("X", "YY", "");
    assertEquals(3, de.size());
    assertEquals(0, de.minLength());
    assertEquals(2, de.maxLength());
    assertEquals("X", de.name(0));
    assertEquals("YY", de.name(1));
    assertEquals("", de.name(2));
  }

  public void test1() {
    final DescriptionCommon de = new DescriptionCommon("YY", "X", "");
    assertEquals(0, de.minLength());
    assertEquals(2, de.maxLength());
  }

  public void test2() {
    final DescriptionCommon de = new DescriptionCommon("", "YY", "X");
    assertEquals(0, de.minLength());
    assertEquals(2, de.maxLength());
  }
}
