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

import static com.rtg.util.StringUtils.LS;

import java.io.IOException;
import java.util.Arrays;

import com.rtg.index.IndexBase.IndexState;
import com.rtg.index.params.CreateParams;
import com.rtg.util.Utils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.Talkback;
import com.rtg.util.test.NanoRegression;

import junit.framework.TestCase;

/**
 */
public class IndexCompressedExtendedTest extends TestCase {

  protected NanoRegression mNano = null;

  @Override
  public void setUp() {
    Diagnostic.setLogStream();
    mNano = new NanoRegression(this.getClass(), false);
  }

  @Override
  public void tearDown() throws Exception {
    // clear the module name so later tests don't report SlimException to the
    // Talkback system
    Talkback.setModuleName(null);
    try {
      mNano.finish();
    } finally {
      mNano = null;
    }
  }

  protected IndexCompressed getIndex(final long size, final int hashBits, final Integer threshold) {
    final CreateParams params = new CreateParams(size, hashBits, hashBits, true, true, false);
    //System.err.println(params);
    return new IndexCompressed(params, threshold, false, threshold, threshold, 2);
  }

  private static long[] extend(final long l) {
    return new long[] {l, 0};
  }

  public final void testStateTwoPass() throws IOException {
    final IndexCompressed hi = getIndex(10L, 67, Integer.MAX_VALUE);
    checkNoSearch(hi);
    assertEquals(IndexState.PRE_ADD, hi.mState);
    hi.add(extend(0), 0);
    hi.freeze();
    assertEquals(IndexState.ADD, hi.mState);
    hi.add(extend(0), 0);
    hi.freeze();
    assertEquals(IndexState.FROZEN, hi.mState);
    //assertEquals(1, hi.searchCount(extend(0)));
  }

  protected void checkNoSearch(final IndexCompressed index) throws IOException {
    try {
      index.search(extend(0), new Finder() {
        @Override
        public void found(final long id) {
          //do nothing
        }
      });
      fail();
    } catch (final IllegalStateException e) {
      //expected
    }
  }

  /**
   * Multiple entries for the same hash.
   * @throws IOException If an I/O error occurs
   */
  public void testScanEmpty() throws IOException {
    final IndexCompressed hi = getIndex(10L, 67, 2);
    try {
      hi.scanAll(null);
      fail();
    } catch (final IllegalStateException e) {

    }
    hi.freeze();
    hi.freeze();
    hi.scanAll(
      new FinderHashValueExtended() {
        @Override
        public void found(long[] hash, long value) {
          //System.err.println("hash=" + hash + " value=" + value);
          fail();
        }
      }
    );
  }

  public final void testTrickySize0() {
    final IndexCompressed hi = getIndex(10L, 67, 2);
    try {
      hi.add(extend(0), 0);
      hi.freeze();
    } catch (final RuntimeException e) {
      assertEquals("Too many items pre-added:1 > 0", e.getMessage());
    }
  }
  public final void testStateSinglePass() throws IOException {
    final IndexCompressed hi = getIndex(10L, 67, Integer.MAX_VALUE);
    checkNoSearch(hi);
    assertEquals(IndexState.PRE_ADD, hi.mState);
    hi.freeze();
    hi.freeze();
    try {
      hi.add(extend(0), 0);
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
    final IndexCompressed hi = getIndex(2L, 65, Integer.MAX_VALUE);
    hi.add(extend(1), 1);
    hi.add(extend(1), 1);
    hi.add(extend(1), 1);
    try {
      hi.freeze();
      fail();
    } catch (final RuntimeException e) {
      assertEquals("Too many items pre-added:3 > 2", e.getMessage());
    }
  }

  public final void testBadConstructor() {
    try {
      getIndex(-1L, 16, Integer.MAX_VALUE);
      fail();
    } catch (final RuntimeException e) {
      //expected
    }
    getIndex(1L, 1, Integer.MAX_VALUE);
    try {
      getIndex(1L, 0, Integer.MAX_VALUE);
      fail();
    } catch (final RuntimeException e) {
      //expected
    }
    try {
      getIndex(1L, 64, 0);
      fail();
    } catch (final RuntimeException e) {
      //expected
    }
  }

  public void testCompressDecompressA1() {
    final IndexCompressed hi = getIndex(2L, 65, Integer.MAX_VALUE);
    checkCompress(hi, new long[] {1, 1});
    checkCompress(hi, new long[] {-1, 0});
  }

  public void testCompressDecompressA2() {
    final IndexCompressed hi = getIndex(2L, 66, Integer.MAX_VALUE);
    checkCompress(hi, new long[] {0, 3});
    checkCompress(hi, new long[] {-1, 3});
  }

  private void checkCompress(IndexCompressed ix, long[] exp) {
    final long u = ix.position(exp);
    final long l = ix.compressHash(exp);
    final long[] actual = ix.decompressHashExtended(u, l);
    if (!Arrays.equals(exp, actual)) {
      final StringBuilder sb = new StringBuilder();
      sb.append(LS);
      sb.append("  exp=").append(Utils.toBits(exp)).append(LS);
      sb.append("  act=").append(Utils.toBits(actual)).append(LS);
      sb.append("lower=").append(Utils.toBits(l, 64)).append(LS);
      sb.append("upper=").append(Utils.toBits(l, 64)).append(LS);
      System.err.println(sb.toString());
      fail(sb.toString());
    }
  }

}
