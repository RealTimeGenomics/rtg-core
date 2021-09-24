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
package com.rtg.metagenomics;

import java.util.Arrays;

import junit.framework.TestCase;

/**
 */
public class SpeciesMapTest extends TestCase {

  public void test() {
    final SpeciesMap sm = new SpeciesMap();
    assertEquals(0, sm.id(1));
    assertEquals(1, sm.id(2));
    sm.put(4, 1);
    assertEquals(0, sm.id(1));
    assertEquals(1, sm.id(2));
    assertEquals(2, sm.id(3));
    assertEquals(1, sm.id(4));
    assertEquals("[1, 2, 3]", Arrays.toString(sm.taxonIds()));
    assertEquals("[1, 2]", Arrays.toString(sm.subset(0, 1).taxonIds()));
  }
}
