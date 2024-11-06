/*
 * Copyright (c) 2018. Real Time Genomics Limited.
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
