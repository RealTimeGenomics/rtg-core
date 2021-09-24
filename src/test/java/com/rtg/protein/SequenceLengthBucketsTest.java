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
