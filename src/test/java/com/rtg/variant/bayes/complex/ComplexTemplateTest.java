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

package com.rtg.variant.bayes.complex;

import java.util.Arrays;

import junit.framework.TestCase;

/**
 */
public class ComplexTemplateTest extends TestCase {

  public void test0() {
    final byte[] template = {};
    final ComplexTemplate cot = new ComplexTemplate(template, "", 0, 0);
    assertEquals(0, cot.getStart());
    assertEquals(0, cot.getEnd());
    assertTrue(template == cot.templateBytes());
    //System.err.println(cot);
    assertEquals(0, cot.getLength());
    assertEquals("", cot.replaceString());
    assertTrue(Arrays.equals(new byte[] {}, cot.replaceBytes()));
  }

  public void test1() {
    final byte[] template = {3, 2, 1, 4};
    final ComplexTemplate cot = new ComplexTemplate(template, "", 3, 4);
    assertEquals(3, cot.getStart());
    assertEquals(4, cot.getEnd());
    assertTrue(template == cot.templateBytes());
    //System.err.println(cot);
    assertEquals(1, cot.getLength());
    assertEquals("T", cot.replaceString());
    assertTrue(Arrays.equals(new byte[] {4}, cot.replaceBytes()));
  }
}
