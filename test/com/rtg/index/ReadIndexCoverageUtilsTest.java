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
package com.rtg.index;

import com.rtg.index.hash.ngs.ReadDecoder;
import com.rtg.index.params.CreateParams;
import com.rtg.util.StringUtils;

import junit.framework.TestCase;

/**
 * test for class
 */
public class ReadIndexCoverageUtilsTest extends TestCase {

  public void testFindRejectedReads() {
    final CreateParams create = new CreateParams(100, 36, 36, false, false, false);
    final IndexSimple index = new IndexSimple(create, new RepeatFrequencyFilterMethod(100, false, 100, 0), 1);
    fillIndex(index,
            new long[] {0, 1, 2, 3, 4, 5, 6, 7, 8, 9},
            new long[] {20, 22, 25, 45, 92, 106, 50, 44, 49, 90});
    //reads 2, 4, 5, 9, 10 exist
    final int[] expected = {0, 1, 3, 6, 7, 8, 11, 12, 13, 14};
    final int[] res = ReadIndexCoverageUtils.findRejectedReads(new Index[] {index}, 15, 10);
    assertEquals(expected.length, res.length);
    for (int i = 0; i < res.length; i++) {
      assertEquals(expected[i], res[i]);
    }
  }

  private static void fillIndex(Index index, long[] hashes, long[] values) {
    for (int i = 0; i < hashes.length; i++) {
      index.add(hashes[i], values[i]);
    }
    index.freeze();
  }

  public void testSummariseRejectedReads() {
    final CreateParams create = new CreateParams(100, 36, 36, false, false, false);
    final IndexSimple index = new IndexSimple(create, new RepeatFrequencyFilterMethod(100, false, 100, 0), 1);
    fillIndex(index,
            new long[] {0, 1, 2, 3, 4, 5, 6, 7, 8, 9},
            new long[] {20, 22, 25, 45, 92, 106, 50, 44, 49, 90});
    final int[] res = ReadIndexCoverageUtils.findRejectedReads(new Index[] {index}, 15, 10);
    assertEquals("Rejected reads total: 10" + StringUtils.LS
            + "Rejected: 0 first: " + true + StringUtils.LS
            + "Rejected: 1 first: " + true + StringUtils.LS
            + "Rejected: 3 first: " + true + StringUtils.LS
            + "Rejected: 6 first: " + true + StringUtils.LS
            + "Rejected: 7 first: " + true + StringUtils.LS
            + "Rejected: 8 first: " + true + StringUtils.LS
            + "Rejected: 11 first: " + true + StringUtils.LS
            + "Rejected: 12 first: " + true + StringUtils.LS
            + "Rejected: 13 first: " + true + StringUtils.LS
            + "Rejected: 14 first: " + true + StringUtils.LS,
            ReadIndexCoverageUtils.summariseRejectedReads(res, ReadDecoder.SINGLE_END));
  }

  public void testSummariseRejectedReadsPaired() {
    final CreateParams create = new CreateParams(100, 36, 36, false, false, false);
    final IndexSimple index = new IndexSimple(create, new RepeatFrequencyFilterMethod(100, false, 100, 0), 1);
    fillIndex(index,
            new long[] {0, 1, 2, 3, 4, 5, 6, 7, 8, 9},
            new long[] {20, 22, 25, 45, 92, 106, 50, 44, 49, 90});
    final int[] res = ReadIndexCoverageUtils.findRejectedReads(new Index[] {index}, 16, 10);
    assertEquals("Rejected arms total: 11" + StringUtils.LS
            + "Rejected: 0 first: " + true + StringUtils.LS
            + "Rejected: 0 first: " + false + StringUtils.LS
            + "Rejected: 1 first: " + false + StringUtils.LS
            + "Rejected: 3 first: " + true + StringUtils.LS
            + "Rejected: 3 first: " + false + StringUtils.LS
            + "Rejected: 4 first: " + true + StringUtils.LS
            + "Rejected: 5 first: " + false + StringUtils.LS
            + "Rejected: 6 first: " + true + StringUtils.LS
            + "Rejected: 6 first: " + false + StringUtils.LS
            + "Rejected: 7 first: " + true + StringUtils.LS
            + "Rejected: 7 first: " + false + StringUtils.LS,
            ReadIndexCoverageUtils.summariseRejectedReads(res, ReadDecoder.PAIRED_END));
  }
}
