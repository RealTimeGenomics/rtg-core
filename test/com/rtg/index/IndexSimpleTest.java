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
