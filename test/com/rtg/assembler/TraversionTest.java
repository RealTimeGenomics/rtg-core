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

package com.rtg.assembler;

import java.util.HashSet;
import java.util.Set;

import junit.framework.TestCase;

/**
 */
public class TraversionTest extends TestCase {
  Set<Long> longs(long... list) {
    Set<Long> result = new HashSet<>();
    for (long l :list) {
      result.add(l);
    }
    return result;
  }

  public void test() {
    Traversion t = new Traversion(longs(1, 2, 3), longs(-1, -2, 3));
    assertEquals(longs(1, 2, 3), t.next());
    assertEquals(longs(-1, -2, 3), t.previous());
  }
}
