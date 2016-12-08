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
import java.util.ArrayList;
import java.util.Arrays;
import java.util.ConcurrentModificationException;
import java.util.List;

import com.rtg.index.IndexBase.IndexState;
import com.rtg.index.params.CreateParams;

/**
 */
public class IndexCompressedTest extends AbstractIndexTest {

  public IndexCompressedTest() {
    super(true);
  }

  @Override
  protected IndexBase getIndex(final long size, final int hashBits, final Integer threshold) {
    return new IndexCompressed(new CreateParams(size, hashBits, hashBits, 31, true, true, false, false), new FixedRepeatFrequencyFilterMethod(threshold), 2);
  }

  public final void testStateTwoPass() throws IOException {
    final IndexBase hi = getIndex(10L, 16, Integer.MAX_VALUE);
    checkNoSearch(hi);
    assertEquals(IndexState.PRE_ADD, hi.mState);
    hi.add(0, 0);
    hi.freeze();
    assertEquals(IndexState.ADD, hi.mState);
    hi.add(0, 0);
    hi.freeze();
    assertEquals(IndexState.FROZEN, hi.mState);
    assertEquals(1, hi.count(0));
  }

  /*
   * Adds the hashes and ids to the index.
   * Checks that they can all be found again and that the missHashes cannot.
   * Does this twice with a melt between.
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
    add(index, hashes, ids);
    index.globalIntegrity();
    //System.err.println("before freeze2: " + index.toString());
    index.freeze();
    //System.err.println("frozen: " + index.toString());
    index.globalIntegrity();
    checkSearch(index, hashes, ids, missHashes, checkCount);
    checkEmpty(index, missHashes, checkCount);

    if (testId != null) {
      mNano.check("indexcompressedtest_" + testId + ".txt", index.toString());
    }
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
    hi.freeze();
    hi.scan(new FinderHashValue() {
      @Override
      public void found(long hash, long value) {
        //System.err.println("hash=" + hash + " value=" + value);
        fail();
      }
    });
  }

  public final void testTrickySize0() {
    final IndexBase hi = getIndex(0, 16, 2);
    try {
      hi.add(0, 0);
      hi.freeze();
    } catch (final RuntimeException e) {
      assertEquals("Too many items pre-added:1 > 0", e.getMessage());
    }
  }
  public final void testStateSinglePass() throws IOException {
    final IndexBase hi = getIndex(10L, 16, Integer.MAX_VALUE);
    checkNoSearch(hi);
    assertEquals(IndexState.PRE_ADD, hi.mState);
    hi.freeze();
    hi.freeze();
    try {
      hi.add(0, 0);
      fail();
    } catch (final IllegalStateException | AssertionError e) {
      //expected
    }
    try {
      hi.freeze();
      fail();
    } catch (final IllegalStateException e) {
      assertEquals("Index closed twice", e.getMessage());
    } catch (final AssertionError e) {
      //expected if assertions turned on
    }
    assertEquals(IndexState.FROZEN, hi.mState);
  }

  public final void testTooMany() {
    final IndexBase hi = getIndex(2L, 16, Integer.MAX_VALUE);
    hi.add(1, 1);
    hi.add(1, 1);
    hi.add(1, 1);
    try {
      hi.freeze();
      fail();
    } catch (final RuntimeException e) {
      assertEquals("Too many items pre-added:3 > 2", e.getMessage());
    }
  }

    static class TestFinder implements FinderHashValue {
    List<Long> mHash = new ArrayList<>();
    List<Long> mValue = new ArrayList<>();

    @Override
    public void found(long hash, long value) {
      mHash.add(hash);
      mValue.add(value);
    }
  }

  public final void testUnsafeError() {
    final CreateParams params = new CreateParams(51, 8, 8, 3, true, true, true, false);
    final IndexExtended countIndex = new IndexCompressed(params, new UnfilteredFilterMethod(), 4);

    final long[] hashes = {48, 192, 3, 48, 192, 3, 48, 192, 1, 4, 19, 77, 53, 104, 88, 96, 61
        , 48, 192, 1, 4, 19, 77, 53, 104, 88, 96, 61, 48, 192, 1, 4, 18, 73, 37, 105, 88, 96, 61, 48
        , 192, 1, 4, 18, 73, 37, 105, 88, 96, 61, 13
    };
    try {
      for (int add = 0; add < 2; ++add) {
        for (long hash : hashes) {
          countIndex.add(hash, 0L);
        }

        countIndex.freeze();
      }
      fail();
    } catch (ConcurrentModificationException e) {
      // Expected
    }
  }
  public final void testAssemblyNonDeterminism() throws IOException {
    for (int i = 0; i < 10; ++i) {
      final CreateParams params = new CreateParams(51, 8, 8, 3, true, true, false, false);
      final IndexExtended countIndex = new IndexCompressed(params, new UnfilteredFilterMethod(), 4);

      final long[] hashes = {48, 192, 3, 48, 192, 3, 48, 192, 1, 4, 19, 77, 53, 104, 88, 96, 61
          , 48, 192, 1, 4, 19, 77, 53, 104, 88, 96, 61, 48, 192, 1, 4, 18, 73, 37, 105, 88, 96, 61, 48
          , 192, 1, 4, 18, 73, 37, 105, 88, 96, 61, 13
      };
      for (int add = 0; add < 2; ++add) {
        for (long hash : hashes) {
          countIndex.add(hash, 0L);
        }
        countIndex.freeze();
      }
      final TestFinder finder = new TestFinder();
      countIndex.scan(finder);
      final List<Long> expected = Arrays.asList(1L, 1L, 1L, 1L, 3L, 3L, 4L, 4L, 4L, 4L, 13L, 18L, 18L, 19L, 19L, 37L, 37L, 48L
          , 48L, 48L, 48L, 48L, 48L, 53L, 53L, 61L, 61L, 61L, 61L, 73L, 73L, 77L, 77L, 88L, 88L, 88L, 88L, 96L, 96L, 96L
          , 96L, 104L, 104L, 105L, 105L, 192L, 192L, 192L, 192L, 192L, 192L);
      assertEquals(expected, finder.mHash);

      assertEquals(73, countIndex.getHash(30));
    }
  }
}
