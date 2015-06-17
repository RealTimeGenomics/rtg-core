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
package com.rtg.ngs.blocking;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;

import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.LogSimple;
import com.rtg.util.io.MemoryPrintStream;

import junit.framework.TestCase;

/**
 * Tests corresponding class
 */
public class ReadBlockerTest extends TestCase {

  public void testBlocking() {
    final ByteArrayOutputStream bos = new ByteArrayOutputStream();
    try (LogSimple ls = new LogSimple(new PrintStream(bos))) {
      Diagnostic.setLogStream(ls);
      for (int k = 1; k < 65536; k <<= 1) {
        final ReadBlocker b = getReadBlocker(2, k);
        for (int j = 0; j < k; j++) {
          assertFalse(b.isBlocked(1));
          assertFalse(b.isBlocked(0));
          b.increment(1);
        }
        assertTrue(b.isBlocked(1));
        assertFalse(b.isBlocked(0));
        b.close();
      }
    } finally {
      Diagnostic.setLogStream();
    }
    final String s = bos.toString();
    assertTrue(s, s.contains("Statistics of "));
    assertTrue(s.contains("Total reads 2"));
    assertTrue(s.contains("1 reads had count 0"));
    assertTrue(s.contains("1 reads had count 64"));
    assertTrue(s.contains(expectedPairingsString()));
    //System.err.println(s);
  }

  protected String expectedPairingsString() {
    return "blocked pairings";
  }

  public void test() {
    final ReadBlocker b = getReadBlocker(1, 255);
    b.increment(0);
    assertEquals(1, b.getCount(0));
    b.reset(0);
    assertEquals(0, b.getCount(0));
    for (int i = 0; i < 65535; i++) {
      assertEquals(i, b.getCount(0));
      b.increment(0);
    }
    assertEquals(65535, b.getCount(0));
    b.increment(0);
    assertEquals(65535, b.getCount(0));
    final MemoryPrintStream ps = new MemoryPrintStream();
    Diagnostic.setLogStream(ps.printStream());
    b.close();
    assertTrue(ps.toString(), ps.toString().contains(">= 65535"));
    Diagnostic.setLogStream();
    ps.close();
  }

  public void testCons() {
    try {
      getReadBlocker(0, 0);
      fail();
    } catch (final IllegalArgumentException e) {
      // ok
    }
    try {
      getReadBlocker(0, 65536);
      fail();
    } catch (final IllegalArgumentException e) {
      // ok
    }
  }

  ReadBlocker getReadBlocker(int numReads, int thres) {
    return new ReadBlocker(numReads, thres);
  }
}
