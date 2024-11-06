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
    final CreateParams create = new CreateParams(100, 36, 36, 31, false, true, false, false);
    final IndexSimple index = new IndexSimple(create, new FixedRepeatFrequencyFilterMethod(100), 1);
    fillIndex(index,
            new long[] {0, 1, 2, 3, 4, 5, 6, 7, 8, 9},
            new long[] {20, 22, 25, 45, 92, 106, 50, 44, 49, 90});
    //reads 2, 4, 5, 9, 10 exist
    final int[] expected = {0, 1, 3, 6, 7, 8, 11, 12, 13, 14};
    final int[] res = ReadIndexCoverageUtils.findRejectedReads(new Index[] {index}, 15, 10);
    assertEquals(expected.length, res.length);
    for (int i = 0; i < res.length; ++i) {
      assertEquals(expected[i], res[i]);
    }
  }

  private static void fillIndex(Index index, long[] hashes, long[] values) {
    for (int i = 0; i < hashes.length; ++i) {
      index.add(hashes[i], values[i]);
    }
    index.freeze();
  }

  public void testSummariseRejectedReads() {
    final CreateParams create = new CreateParams(100, 36, 36, 31, false, true, false, false);
    final IndexSimple index = new IndexSimple(create, new FixedRepeatFrequencyFilterMethod(100), 1);
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
    final CreateParams create = new CreateParams(100, 36, 36, 31, false, true, false, false);
    final IndexSimple index = new IndexSimple(create, new FixedRepeatFrequencyFilterMethod(100), 1);
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
