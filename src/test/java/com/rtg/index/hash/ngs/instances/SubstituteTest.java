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
package com.rtg.index.hash.ngs.instances;




import java.io.IOException;
import java.io.StringWriter;

import com.rtg.index.hash.ngs.HashFunctionFactory;
import com.rtg.index.hash.ngs.NgsHashFunction;
import com.rtg.index.hash.ngs.ReadCall;
import com.rtg.index.hash.ngs.TemplateCall;
import com.rtg.launcher.HashingRegion;
import com.rtg.util.StringUtils;

public final class SubstituteTest extends AbstractSplitTest {

  @Override
  protected NgsHashFunction getHashFunction(final ReadCall readCall, final TemplateCall templateCall) {
    return null;
  }

  static final class HashFunctionMock implements NgsHashFunction, Cloneable {
    private final Appendable mOut;
    final StringBuilder mSB = new StringBuilder();

    HashFunctionMock(final Appendable out) {
      mOut = out;
    }

    @Override
    public void readAll(final int readId, final boolean reverse) {
    }

    @Override
    public int readLength() {
      return 0;
    }

    @Override
    public void setReadSequences(final long numberReads) {
    }

    @Override
    public void setValues(final int id2, final boolean reverse) {
    }

    @Override
    public int fastScore(final int readId) {
      return 0;
    }

    @Override
    public int indelScore(final int readId) {
      return 0;
    }

    @Override
    public void templateForward(final int endPosition) {
      append(mOut, mSB.toString() + StringUtils.LS);
    }

    @Override
    public void templateBidirectional(final int endPosition) {
      append(mOut, mSB.toString() + StringUtils.LS);
    }

    @Override
    public void templateReverse(final int endPosition) {
      throw new UnsupportedOperationException();
    }

    @Override
    public void endSequence() {
      // do nothing
    }

    @Override
    public void templateSet(final long name, final int length) {
    }

    @Override
    public void hashStep(final byte code) {
      mSB.append(CHARS[code + 1]);
    }

    @Override
    public void hashStep() {
      hashStep((byte) 0);
    }

    @Override
    public int numberWindows() {
      return 0;
    }

    @Override
    public void reset() {
      mSB.delete(0, mSB.length());
    }

    @Override
    public int windowSize() {
      return 0;
    }

    @Override
    public NgsHashFunction threadClone(final HashingRegion region) {
      if (region != HashingRegion.NONE) {
        throw new UnsupportedOperationException();
      }
      try {
        return (NgsHashFunction) super.clone();
      } catch (final CloneNotSupportedException e) {
        throw new RuntimeException(e);
      }
    }

    @Override
    public void threadFinish() {
    }

    /**
     * Standard semantics of clone.
     */
    @Override
    public HashFunctionMock clone() throws CloneNotSupportedException {
      return (HashFunctionMock) super.clone();
    }

    @Override
    public void logStatistics() {
      // do nothing
    }
  }

  /**
   * Test the substitution testing code.
   */
  private void checkSubstitute(final int n, final int expected, final String expStr) throws IOException {
    final StringWriter sb = new StringWriter();
    final HashFunctionFactory hff = new HashFunctionFactory() {
      @Override
      public NgsHashFunction create(final ReadCall readCall, final TemplateCall templateCall) {
        return new HashFunctionMock(sb);
      }

      @Override
      public int hashBits() {
        return 0;
      }

      @Override
      public int numberWindows() {
        return 0;
      }

      @Override
      public int windowBits() {
        return 0;
      }

      @Override
      public int windowSize() {
        throw new UnsupportedOperationException("Not supported yet.");
      }

    };
    final Substitute sub = new Substitute("acg", hff, false);
    sub.substituteProtected(n);
    if (expStr != null) {
      assertEquals(expStr, sb.toString());
    }
    final int actualCount = sb.toString().split(StringUtils.LS).length;
    assertEquals(expected, actualCount);
  }

  private static final String EXPECTED_1 = ""
    + "aag" + StringUtils.LS
    + "aca" + StringUtils.LS
    + "acg" + StringUtils.LS
    + "ccg" + StringUtils.LS;

  public void testSubstitute() throws IOException {
    checkSubstitute(0, 1, "acg" + StringUtils.LS);
    checkSubstitute(1, 4, EXPECTED_1);
    checkSubstitute(2, 7, null);
  }
}

