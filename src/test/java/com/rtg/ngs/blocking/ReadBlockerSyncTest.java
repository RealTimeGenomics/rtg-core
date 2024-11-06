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

import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.MemoryPrintStream;

/**
 * Test class
 */
public class ReadBlockerSyncTest extends ReadBlockerTest {

  @Override
  ReadBlocker getReadBlocker(int reads, int threshold) {
    return new ReadBlockerSync(reads, threshold, "test-sync");
  }

  @Override
  protected String expectedPairingsString() {
    return "test-sync";
  }

  @Override
  public void test() {
    final ReadBlocker b = getReadBlocker(1, 255);
    b.increment(0);
    assertEquals(1, b.getCount(0));
    b.reset(0);
    assertEquals(0, b.getCount(0));
    for (int i = 0; i < 65535; ++i) {
      b.increment(0);
    }
    final MemoryPrintStream ps = new MemoryPrintStream();
    Diagnostic.setLogStream(ps.printStream());
    b.close();
    assertTrue(ps.toString(), ps.toString().contains("1 reads had count 255"));
    Diagnostic.setLogStream();
    ps.close();
  }
}
