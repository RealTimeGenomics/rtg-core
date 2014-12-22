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
package com.rtg.variant.bayes.multisample.lineage;

import junit.framework.TestCase;

/**
 * Test class.
 */
public class VariableTest extends TestCase {

  public void test() {
    final Variable v = new Variable("X", 42);
    assertEquals("X", v.toString());
    assertEquals(42, v.size());
    assertEquals(v, v);
    assertEquals(v, new Variable("X", 42));
  }
}
