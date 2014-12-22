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
package com.rtg.util;

import junit.framework.TestCase;


/**
 * Tests corresponding class
 */
public class ComparablePairTest extends TestCase {

  public void testCompare() {
    final ComparablePair<String, String> p1 = new ComparablePair<>("a", "b");
    assertEquals(0, p1.compareTo(new ComparablePair<>("a", "b")));
    assertEquals(-1, p1.compareTo(new ComparablePair<>("b", "b")));
    assertEquals(1, p1.compareTo(new ComparablePair<>("a", "a")));
  }
}
