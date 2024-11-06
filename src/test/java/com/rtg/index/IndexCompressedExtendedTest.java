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
    final CreateParams params = new CreateParams(size, hashBits, hashBits, 31, true, true, true, false);
    //System.err.println(params);
    return new IndexCompressed(params, new FixedRepeatFrequencyFilterMethod(threshold), 2);
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
        public boolean found(final long id) {
          return true;
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
  public void testScanEmpty() {
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
