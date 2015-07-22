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

import java.io.IOException;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashSet;
import java.util.Set;

import com.rtg.mode.SequenceType;
import com.rtg.reader.MockSequencesReader;

import junit.framework.TestCase;

/**
 *         Date: 30/09/11
 *         Time: 9:54 AM
 */
public class SequenceLengthBucketsTest extends TestCase {
  private static class MockLengthsReader extends MockSequencesReader {
    private final int[] mLengths;
    public MockLengthsReader(int[] lengths) {
      super(SequenceType.PROTEIN, lengths.length);
      mLengths = lengths;
    }
    @Override
    public int[] sequenceLengths(long start, long end) {
      return Arrays.copyOfRange(mLengths, (int) start, (int) end);
    }

  }

  public void testBucketing() throws Exception {
    final int[] lengths = {100, 101, 102, 103, 104, 23, 105, 106, 107, 50, 108};
    final SequenceLengthBuckets buckets = new SequenceLengthBuckets(new MockLengthsReader(lengths), 99);
    assertEquals(4, buckets.getBuckets().size());
    assertEquals(2, buckets.getCount(99));
    assertEquals(3, buckets.getCount(102));
    assertEquals(3, buckets.getCount(105));
    assertEquals(1, buckets.getCount(108));
    assertEquals(0, buckets.getCount(9999));
  }

  public void testBucketList() throws IOException {
    final int[] lengths = {100, 101, 102, 103, 104, 105, 24, 106, 107, 108};
    final SequenceLengthBuckets buckets = new SequenceLengthBuckets(new MockLengthsReader(lengths), 99);
    final Collection<Long> b = buckets.getBuckets();
    final Set<Long> expected = new HashSet<>();
    expected.addAll(Arrays.asList(99L, 102L, 105L, 108L));
    assertEquals(expected, b);
  }
}
