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

package com.rtg.variant.sv;

import com.rtg.util.TestUtils;
import com.rtg.util.io.MemoryPrintStream;

import junit.framework.TestCase;

/**
 */
public class SignalSumTest extends TestCase {

  public void test() {
    final SamCounts sa = new SamArray(10);
    sa.increment(0);
    sa.increment(4);
    sa.increment(4);
    final Signal sigc = new SignalCount(sa, -3, 3, "");
    try (MemoryPrintStream ps = new MemoryPrintStream()) {
      final Signal sig = new SignalSum(ps.printStream(), "blah", sigc, sigc);
      assertEquals("Sum:", sig.toString());
      assertEquals(1.0 * 2, sig.value(0));
      assertEquals(1.0 * 2, sig.value(1));
      assertEquals(3.0 * 2, sig.value(2));
      assertEquals(3.0 * 2, sig.value(3));
      assertEquals(2.0 * 2, sig.value(4));
      assertEquals(2.0 * 2, sig.value(5));
      assertEquals(2.0 * 2, sig.value(6));
      assertEquals(2.0 * 2, sig.value(7));
      assertEquals(0.0 * 2, sig.value(8));
      assertEquals(0.0 * 2, sig.value(9));
      assertEquals(0.0 * 2, sig.value(10));
      assertEquals("blah", sig.columnLabel());
      TestUtils.containsAll(ps.toString().replace("\t", "  "), "0  1.00  1.00  2.00", "10  0.00  0.00  0.00");
    }
  }

  public void test0() {
    final Signal sig = new SignalSum("blah");
    assertEquals(0.0, sig.value(0));
    assertEquals("blah", sig.columnLabel());
  }

}
