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
package com.rtg.index.similarity;

import com.rtg.index.FixedRepeatFrequencyFilterMethod;
import com.rtg.index.Index;
import com.rtg.index.params.CreateParams;
import com.rtg.util.StringUtils;
import com.rtg.util.diagnostic.Diagnostic;

import junit.framework.TestCase;

/**
 */
public class IndexSimilarityTest extends TestCase {

  @Override
  public void setUp() {
    Diagnostic.setLogStream();
  }

  /**
   * Adds the hashes and ids to the index.
   * Checks that they can all be found again and that the missHashes cannot.
   * @throws IllegalStateException
   */
  private void check(final IndexSimilarity index, final long[] hashes, final int[] ids, final String expected) {
    index.globalIntegrity();
    assertEquals(hashes.length, ids.length);
    add(index, hashes, ids);
    index.globalIntegrity();
    index.freeze();
    add(index, hashes, ids);
    index.globalIntegrity();
    index.freeze();
    index.globalIntegrity();
    //System.err.println(index.toString());
    assertEquals(expected, checkSimilarity(index));
  }

  /**
   * Add the hashes and ids to the index.
   */
  private void add(final Index index, final long[] hashes, final int[] ids) {
    for (int i = 0; i < hashes.length; ++i) {
      //System.err.println(hashes[i] + ":" + ids[i]);
      index.add(hashes[i], ids[i]);
    }
  }

  /**
   * Check that the hashes and ids can be found in the index.
   */
  private String checkSimilarity(final IndexSimilarity index) {
    return index.similarity(10).toString();
  }

  public final void testSimilarity1() {
    final CreateParams params = new CreateParams(100, 32, 32, 31, true, true, false, false);
    final IndexSimilarity index = new IndexSimilarity(params, new FixedRepeatFrequencyFilterMethod(6), false, 1);
    final String expected = ""
      + "SimilarityMatrix 10" + StringUtils.LS
      + "[0]\t1\t1\t3\t1\t0\t0\t0\t0\t0\t0" + StringUtils.LS
      + "[1]\t1\t1\t3\t1\t0\t0\t0\t0\t0\t0" + StringUtils.LS
      + "[2]\t3\t3\t9\t3\t0\t0\t0\t0\t0\t0" + StringUtils.LS
      + "[3]\t1\t1\t3\t5\t6\t0\t0\t0\t0\t0" + StringUtils.LS
      + "[4]\t0\t0\t0\t6\t9\t0\t0\t0\t0\t0" + StringUtils.LS
      + "[5]\t0\t0\t0\t0\t0\t1\t0\t0\t0\t0" + StringUtils.LS
      + "[6]\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0" + StringUtils.LS
      + "[7]\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0" + StringUtils.LS
      + "[8]\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0" + StringUtils.LS
      + "[9]\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0" + StringUtils.LS
      + "" + StringUtils.LS
      ;
    check(index,
          new long[] {1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3},
          new int[]  {1, 2, 3, 0, 2, 2, 3, 4, 3, 4, 4, 5},
          expected);
  }

  public final void testSimilarity2() {
    final CreateParams params = new CreateParams(100, 32, 32, 31, true, true, false, false);
    final IndexSimilarity index = new IndexSimilarity(params, new FixedRepeatFrequencyFilterMethod(6), false, 1);
    final String expected = ""
      + "SimilarityMatrix 10" + StringUtils.LS
      + "[0]\t1\t1\t3\t1\t0\t0\t0\t0\t0\t0" + StringUtils.LS
      + "[1]\t1\t1\t3\t1\t0\t0\t0\t0\t0\t0" + StringUtils.LS
      + "[2]\t3\t3\t9\t3\t0\t0\t0\t0\t0\t0" + StringUtils.LS
      + "[3]\t1\t1\t3\t1\t0\t0\t0\t0\t0\t0" + StringUtils.LS
      + "[4]\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0" + StringUtils.LS
      + "[5]\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0" + StringUtils.LS
      + "[6]\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0" + StringUtils.LS
      + "[7]\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0" + StringUtils.LS
      + "[8]\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0" + StringUtils.LS
      + "[9]\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0" + StringUtils.LS
      + "" + StringUtils.LS
      ;
    check(index,
          new long[] {1, 1, 1, 1, 1, 1, },
          new int[]  {1, 2, 3, 0, 2, 2, },
          expected);
  }

}
