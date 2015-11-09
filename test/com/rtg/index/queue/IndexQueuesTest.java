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
package com.rtg.index.queue;

import static com.rtg.util.StringUtils.LS;

import java.io.IOException;
import java.io.PrintStream;

import com.rtg.index.Add;
import com.rtg.index.Finder;
import com.rtg.index.FinderHashValue;
import com.rtg.index.Index;
import com.rtg.util.diagnostic.Diagnostic;

import junit.framework.TestCase;


/**
 */
public class IndexQueuesTest extends TestCase {

  private static class MockIndex implements Index {

    private static final long MASK = (1L << 3) - 1;

    private final StringBuilder mSb = new StringBuilder();

    @Override
    public synchronized void add(final long hash, final long id) {
      mSb.append("add radix=").append(hash >> 3).append(" hash=").append(hash & MASK).append(" id=").append(id).append(LS);
    }

    @Override
    public long bytes() {
      return 0;
    }

    @Override
    public void dumpValues(final PrintStream out) {
    }

    @Override
    public synchronized void freeze() {
      mSb.append("freeze").append(LS);
    }

    @Override
    public long getHash(final long found) {
      return 0;
    }

    @Override
    public long getValue(long found) {
      return 0;
    }


    @Override
    public String infoString() {
      return null;
    }

    @Override
    public long numberEntries() {
      return 0;
    }

    @Override
    public long numberHashes() {
      return 0;
    }

    @Override
    public String perfString() {
      return null;
    }

    @Override
    public void search(final long hash, final Finder finder) throws IOException, IllegalStateException {
    }

    @Override
    public void scan(FinderHashValue finder) throws IOException, IllegalStateException {
    }

    @Override
    public long contains(final long hash) {
      return 0;
    }

    @Override
    public int count(final long hash) {
      return 0;
    }

    @Override
    public String toString() {
      return mSb.toString();
    }

    @Override
    public int maxHashCount() {
      throw new UnsupportedOperationException("Not supported yet.");
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
    final Index ix = new MockIndex();
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
      assertEquals("0", e.getMessage());
    }
  }
}
