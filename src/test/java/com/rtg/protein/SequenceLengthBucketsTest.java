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
package com.rtg.protein;

import java.util.Arrays;
import java.util.Collection;
import java.util.HashSet;
import java.util.Set;

import com.rtg.AbstractTest;

/**
 */
public class SequenceLengthBucketsTest extends AbstractTest {

  public void testBucketing() {
    final int[] lengths = {100, 101, 102, 103, 104, 23, 105, 106, 107, 50, 108};
    final SequenceLengthBuckets buckets = new SequenceLengthBuckets(3, 99, lengths);
    assertEquals(4, buckets.getBuckets().size());
    assertEquals(2, buckets.getCount(99));
    assertEquals(3, buckets.getCount(102));
    assertEquals(3, buckets.getCount(105));
    assertEquals(1, buckets.getCount(108));
    assertEquals(0, buckets.getCount(9999));
  }

  public void testBucketList() {
    final int[] lengths = {100, 101, 102, 103, 104, 105, 24, 106, 107, 108};
    final SequenceLengthBuckets buckets = new SequenceLengthBuckets(3, 99, lengths);
    final Collection<Long> b = buckets.getBuckets();
    final Set<Long> expected = new HashSet<>(Arrays.asList(99L, 102L, 105L, 108L));
    assertEquals(expected, b);
  }
}
