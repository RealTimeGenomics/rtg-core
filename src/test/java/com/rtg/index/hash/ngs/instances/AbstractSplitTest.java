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
import java.util.Collection;

import com.rtg.index.hash.ngs.NgsHashFunction;
import com.rtg.index.hash.ngs.ReadCall;
import com.rtg.index.hash.ngs.TemplateCall;
import com.rtg.launcher.AbstractNanoTest;
import com.rtg.launcher.HashingRegion;
import com.rtg.mode.DNA;
import com.rtg.util.MultiMap;
import com.rtg.util.Pair;
import com.rtg.util.StringUtils;
import com.rtg.util.TestUtils;
import com.rtg.util.Utils;
import com.rtg.util.integrity.Exam;



/**
 */
public abstract class AbstractSplitTest extends AbstractNanoTest {

  static void append(final Appendable ap, final String str) {
    try {
      ap.append(str);
    } catch (final IOException e) {
      fail();
    }
  }

  /*
    ((c & 15) >>> 1) ^ 3
    gives g->0, t->1, c->2, a->3, n->4
   */

  private static final byte[] ORDINALS = {
    (byte) (DNA.G.ordinal() - 1),
    (byte) (DNA.T.ordinal() - 1),
    (byte) (DNA.C.ordinal() - 1),
    (byte) (DNA.A.ordinal() - 1),
    (byte) (DNA.N.ordinal() - 1),
  };

  public static void encode(final NgsHashFunction hf, final String dna) {
    for (int i = 0; i < dna.length(); ++i) {
      encode(hf, dna, i);
    }
  }

  public static void encode(final NgsHashFunction hf, final String dna, int pos) {
    hf.hashStep(ORDINALS[((dna.charAt(pos) & 15) >>> 1) ^ 3]);
  }

  /**
   * Mock <code>ReadCall</code> and accumulate results.
   */
  public static class ReadCallMock implements ReadCall {
    private final Appendable mOut;

    /**
     * @param out output
     */
    public ReadCallMock(final Appendable out) {
      mOut = out;
    }
    @Override
    public void readCall(final int id, final long hash, final int index) {
      append(mOut, "readCall index=" + index + " id=" + id + StringUtils.LS);
      append(mOut, " hash=" + Utils.toBitsSep(hash) + StringUtils.LS);
    }
  }


  static class ReadCallAccumulate implements ReadCall {
    private final MultiMap<Pair<Long, Integer>, Integer> mMap = new MultiMap<>();
    @Override
    public void readCall(final int id, final long hash, final int index) {
      mMap.put(new Pair<>(hash, index), id);
    }

    MultiMap<Pair<Long, Integer>, Integer> map() {
      return mMap;
    }

  }

  /**
   * Mock <code>TemplateCall</code> and accumulate results.
   */
  public static class TemplateCallMock implements TemplateCall, Cloneable {
    private final Appendable mOut;
    /**
     * @param out output
     */
    public TemplateCallMock(final Appendable out) {
      mOut = out;
    }
    @Override
    public void done() {
      append(mOut, "done" + StringUtils.LS);
    }
    @Override
    public void endSequence() {
      append(mOut, "endSequence" + StringUtils.LS);
    }
    @Override
    public void set(final long name, final int length) {
      append(mOut, "set name=" + name + " length=" + length + StringUtils.LS);
    }

    @Override
    public void setReverse(final boolean reverse) {
      assertEquals(false, reverse);
    }

    @Override
    public boolean isReverse() {
      return false;
    }
    @Override
    public void setHashFunction(final NgsHashFunction hashFunction) {
      //do nothing
    }
    @Override
    public void templateCall(final int endPosition, final long hash, final int index) {
      try {
        mOut.append("templateCall position=");
        mOut.append(String.valueOf(endPosition));
        mOut.append(" index=");
        mOut.append(String.valueOf(index));
        mOut.append(StringUtils.LS);
        mOut.append(" hash=");
        mOut.append(String.valueOf(Utils.toBitsSep(hash)));
        mOut.append(StringUtils.LS);
      } catch (final IOException e) {
        fail();
      }
    }
    @Override
    public TemplateCall threadClone(final HashingRegion region) {
      if (region != HashingRegion.NONE) {
        throw new UnsupportedOperationException();
      }
      try {
        return clone();
      } catch (final CloneNotSupportedException e) {
        throw new RuntimeException(e);
      }
    }

    @Override
    public void threadFinish() {
    }

    /**
     */
    @Override
    public TemplateCallMock clone() throws CloneNotSupportedException {
      return (TemplateCallMock) super.clone();
    }
    @Override
    public void logStatistics() {
      // do nothing
    }
  }

  static class TemplateCallCheck implements TemplateCall, Cloneable {
    private final MultiMap<Pair<Long, Integer>, Integer> mMap;

    boolean mFound = false;

    boolean mFoundDone = false;

    TemplateCallCheck(final MultiMap<Pair<Long, Integer>, Integer> map) {
      mMap = map;
    }

    @Override
    public boolean isReverse() {
      return false;
    }

    @Override
    public void done() {
      mFoundDone = mFound;
    }
    @Override
    public void endSequence() {
    }
    @Override
    public void set(final long name, final int length) {
    }

    @Override
    public void setReverse(final boolean reverse) {
      // do nothing
    }

    @Override
    public void setHashFunction(final NgsHashFunction hashFunction) {
    }
    @Override
    public void templateCall(final int endPosition, final long hash, final int index) {
      final Collection<Integer> coll = mMap.get(new Pair<>(hash, index));
      mFound |= coll != null && !coll.isEmpty();
    }
    @Override
    public TemplateCall threadClone(final HashingRegion region) {
      if (region != HashingRegion.NONE) {
        throw new UnsupportedOperationException();
      }
      try {
        return clone();
      } catch (final CloneNotSupportedException e) {
        throw new RuntimeException(e);
      }
    }

    @Override
    public void threadFinish() {
    }
    /**
     */
    @Override
    public TemplateCallCheck clone() throws CloneNotSupportedException {
      return (TemplateCallCheck) super.clone();
    }
    @Override
    public void logStatistics() {
      // do nothing
    }
  }

  /**
   * @return a hash function with the hash call set to the identity.
   */
  protected abstract NgsHashFunction getHashFunction(final ReadCall readCall, final TemplateCall templateCall);

  protected void checkr(final String dnar, final String dnat, final String expected) throws IOException {
    final StringWriter sb = new StringWriter();
    final ReadCall rcall = new ReadCallMock(sb);
    final TemplateCall call = new TemplateCallMock(sb);
    final NgsHashFunction hf = getHashFunction(rcall, call);
    Exam.integrity(hf);
    sb.append(hf.toString())
      .append(StringUtils.LS)
      .append("number windows=")
      .append(String.valueOf(hf.numberWindows()))
      .append(StringUtils.LS);
    hf.templateSet(1234, 23);
    encode(hf, dnar);
    hf.readAll(5, false);
    encode(hf, dnat);
    hf.templateForward(7);
    TestUtils.assertResourceEquals("com/rtg/index/hash/ngs/instances/resources/" + expected, sb.toString());
  }

  protected void check(final String dnar, final String dnat, final String expected) throws IOException {
    final StringWriter sb = new StringWriter();
    final ReadCall rcall = new ReadCallMock(sb);
    final TemplateCall call = new TemplateCallMock(sb);
    final NgsHashFunction hf = getHashFunction(rcall, call);
    Exam.integrity(hf);
    sb.append(hf.toString())
      .append(StringUtils.LS)
      .append("number windows=")
      .append(String.valueOf(hf.numberWindows()))
      .append(StringUtils.LS);
    hf.templateSet(1234, 23);
    encode(hf, dnar);
    hf.readAll(5, false);
    encode(hf, dnat);
    hf.templateForward(7);
    assertEquals(expected, sb.toString());
  }

  protected void checkNano(final String id, final String dnar, final String dnat) throws IOException {
    final StringWriter sb = new StringWriter();
    final ReadCall rcall = new ReadCallMock(sb);
    final TemplateCall call = new TemplateCallMock(sb);
    final NgsHashFunction hf = getHashFunction(rcall, call);
    Exam.integrity(hf);
    sb.append(hf.toString())
      .append(StringUtils.LS)
      .append("number windows=")
      .append(String.valueOf(hf.numberWindows()))
      .append(StringUtils.LS);
    hf.templateSet(1234, 23);
    encode(hf, dnar);
    hf.readAll(5, false);
    encode(hf, dnat);
    hf.templateForward(7);
    mNano.check(id, sb.toString());
  }

  static final char[] CHARS = {'n', 'a', 'c', 'g', 't'};
}

