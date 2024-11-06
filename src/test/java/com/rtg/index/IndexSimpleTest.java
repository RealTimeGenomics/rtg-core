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

import java.io.IOException;

import com.rtg.index.params.CreateParams;

/**
 */
public class IndexSimpleTest extends AbstractIndexTest {

  public IndexSimpleTest() {
    super(false);
  }

  /*
   * Adds the hashes and ids to the index.
   * Checks that they can all be found again and that the missHashes cannot.
   * @param expected what toString should show after the build
   */
  @Override
  protected void check(String testId, final IndexBase index, final long[] hashes, final int[] ids, final long[] missHashes, final boolean checkCount) throws IOException {
    index.globalIntegrity();
    //System.out.println("check1 with twopass=" + twoPass);
    assertEquals(hashes.length, ids.length);
    add(index, hashes, ids);
    index.globalIntegrity();
    //System.err.println("before freeze1: " + index.toString());
    index.freeze();
    //System.err.println("frozen: " + index.toString());
    index.globalIntegrity();
    checkSearch(index, hashes, ids, missHashes, checkCount);
    checkEmpty(index, missHashes, checkCount);
  }

  @Override
  protected IndexBase getIndex(final long size, final int hashBits, final Integer threshold) {
    return new IndexSimple(new CreateParams(size, hashBits, hashBits, 31, false, true, false, false), new FixedRepeatFrequencyFilterMethod(threshold), 1);
  }

  /**
   * Multiple entries for the same hash.
   * @throws IOException If an I/O error occurs
   */
  public void testScanEmpty() throws IOException {
    final IndexBase hi = getIndex(10L, 16, 2);
    try {
      hi.scan(null);
      fail();
    } catch (final IllegalStateException e) {

    }
    hi.freeze();
    hi.scan(new FinderHashValue() {
      @Override
      public void found(long hash, long value) {
        System.err.println("hash=" + hash + " value=" + value);
        fail();
      }
    });
  }

  public final void testTrickySize0() {
    final IndexBase hi = getIndex(0, 16, 2);
    try {
      hi.add(0, 0);
      fail();
    } catch (final RuntimeException e) {
      assertEquals("Too many items added:0", e.getMessage());
    }
    hi.freeze();
  }

  private static final int[] EXPECTED_MULTIPLIER = {
    1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, //16
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, //32
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, //64
  };

  public final void testHashSize() {
    for (int i = 1; i <= 64; ++i) {
      final IndexBase hi = getIndex(1, i, 2);
      assertTrue(i + hi.infoString(), hi.infoString().contains(EXPECTED_MULTIPLIER[i - 1] + "\t1\tHash"));
    }
  }

}
