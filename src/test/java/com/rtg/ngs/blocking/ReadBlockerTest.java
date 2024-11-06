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
        for (int j = 0; j < k; ++j) {
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
    for (int i = 0; i < 65535; ++i) {
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
