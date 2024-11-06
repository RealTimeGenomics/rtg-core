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
package com.rtg.index.queue;

import static com.rtg.util.StringUtils.LS;

import com.rtg.index.Add;
import com.rtg.index.Index;
import com.rtg.index.IndexCompressed;
import com.rtg.index.UnfilteredFilterMethod;
import com.rtg.index.params.CreateParams;
import com.rtg.position.MockIndex;
import com.rtg.util.diagnostic.Diagnostic;

import junit.framework.TestCase;

/**
 */
public class IndexQueuesTest extends TestCase {

  private static class MyMockIndex extends MockIndex {

    private static final long MASK = (1L << 3) - 1;

    private final StringBuilder mSb = new StringBuilder();

    @Override
    public synchronized void add(final long hash, final long id) {
      mSb.append("add radix=").append(hash >> 3).append(" hash=").append(hash & MASK).append(" id=").append(id).append(LS);
    }

    @Override
    public synchronized void freeze() {
      mSb.append("freeze").append(LS);
    }

    @Override
    public String toString() {
      return mSb.toString();
    }
  }

  public void test() {
    Diagnostic.setLogStream();
    final IndexQueues iq = new IndexQueues(4, 13, 10, 2, 9);
    assertEquals("IndexQueues: threads=4 radixBits=9 radixSize=3 lowerBits=4", iq.toString());
    iq.integrity();
    final Add a0 = iq.queue(0);
    final Add a1 = iq.queue(1);
    //final Add a2 = iq.queue(2);
    final Add a3 = iq.queue(3);
    final int h1 = (1 << 3) + 1;
    a0.add(h1, 1);
    a1.add(h1, 2);
    final int h2 = (2 << 3) + 2;
    a1.add(h2, 1);
    a0.add(h2, 2);
    a3.add(h2, 3);
    final Index ix = new MyMockIndex();
    iq.freeze(ix);
    final String exp = ""
        + "add radix=1 hash=1 id=1" + LS
        + "add radix=1 hash=1 id=2" + LS
        + "add radix=2 hash=2 id=2" + LS
        + "add radix=2 hash=2 id=1" + LS
        + "add radix=2 hash=2 id=3" + LS
        + "freeze" + LS
        + "add radix=1 hash=1 id=1" + LS
        + "add radix=1 hash=1 id=2" + LS
        + "add radix=2 hash=2 id=2" + LS
        + "add radix=2 hash=2 id=1" + LS
        + "add radix=2 hash=2 id=3" + LS
        + "freeze" + LS
        ;
    assertEquals(exp, ix.toString());
  }

  public void testValues() {
    Diagnostic.setLogStream();
    check(4, 13, 10, 3, 10, 3);
    check(4, 13, 11, 3, 10, 3);
    check(4, 13, 12, 3, 10, 3);
    check(4, 13, 13, 3, 10, 4);
    check(4, 11, 13, 1, 10, 4);
    check(4, 10, 13, 0, 10, 4);
    check(4,  9, 13, 0,  9, 4);
  }

  private void check(final int threads, final int hashBits, final long size, final int xLowerBits, final int xRadixBits, final long xRadixSize) {
    final IndexQueues iq = new IndexQueues(threads, hashBits, size, 0, hashBits);
    iq.integrity();
    assertEquals("IndexQueues: threads=" + threads + " radixBits=" + xRadixBits + " radixSize=" + xRadixSize + " lowerBits=" + xLowerBits, iq.toString());
  }

  public void testBadThreads() {
    try {
      new IndexQueues(0, 13, 10, 0, 2);
    } catch (final RuntimeException e) {
      assertEquals("threads=0", e.getMessage());
    }
  }

  public void testLongs() {
    final CreateParams indexParams = new CreateParams(1, 20, 20, 33, true, true, false, false);
    final IndexQueues iq = new IndexQueues(4, indexParams.hashBits(), indexParams.size(), indexParams.valueBits(), indexParams.initialPointerBits());
    iq.queue(0).add(2, 0b111111111111111111111111111111111L);
    final Index i = new IndexCompressed(indexParams, new UnfilteredFilterMethod(), 1);
    iq.freeze(i);
    assertEquals(0b111111111111111111111111111111111L, i.getValue(i.first(2)));
  }
}
