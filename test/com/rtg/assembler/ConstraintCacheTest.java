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

import java.util.Arrays;

import junit.framework.TestCase;

/**
 */
public class ConstraintCacheTest extends TestCase {
  public void testCache() {
    final ConstraintCache cache = new ConstraintCache();
    cache.addConstraint(4, 5, 199, 2, 200);
    cache.addConstraint(4, 5, 101, 100, 200);
    cache.addConstraint(5, 9, 50, 50, 250);
    cache.addConstraint(5, 9, 40, 60, 200);
    cache.addConstraint(9, 5, 45, 55, 333);

    assertEquals(Arrays.asList(150, 100, 233), cache.find(5, 9).mConstraint);

    assertEquals(Arrays.asList(150, 100, 233), cache.find(9, 5).mConstraint);
    assertNull(cache.find(400));
    assertNull(cache.find(400, 600));

    final ConstraintCache cache2 = new ConstraintCache();
    cache.addConstraint(5, 9, 52, 52, 200);
    cache.addConstraint(5, 9, 42, 62, 204);
    cache.addConstraint(6, 8, 42, 62, 150);

    final ConstraintCache combined = ConstraintCache.combineCaches(Arrays.asList(cache, cache2));
    assertEquals(Arrays.asList(150, 100, 233, 96, 100), combined.find(5, 9).mConstraint);

    assertEquals(Arrays.asList(-1, -1), combined.find(4, 5).mConstraint);
    assertEquals(Arrays.asList(46), combined.find(6, 8).mConstraint);
  }
}
