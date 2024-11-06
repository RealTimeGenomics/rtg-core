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

import static com.rtg.util.StringUtils.LS;

import com.rtg.util.TestUtils;
import com.rtg.util.io.MemoryPrintStream;

import junit.framework.TestCase;

/**
 */
public class AllCountsTest extends TestCase {

  public void test() {
    final AllCounts ac = new AllCounts(5);
    ac.globalIntegrity();
    assertEquals("AllCounts:5", ac.toString());
    TestUtils.containsAll(plot(ac), "0 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000", "4 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000");

    ac.properLeft().increment(0);
    ac.discordantLeft().increment(1, 2);
    ac.unmatedLeft().increment(2, 3);
    ac.properRight().increment(3);
    ac.discordantRight().increment(4, 2);
    ac.unmatedRight().increment(3, 3);
    ac.unpaired().increment(2, 4);
    final String exp = ""
      + "0 1.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000" + LS
      + "1 0.0000 2.0000 0.0000 0.0000 0.0000 0.0000 0.0000" + LS
      + "2 0.0000 0.0000 3.0000 0.0000 0.0000 0.0000 4.0000" + LS
      + "3 0.0000 0.0000 0.0000 1.0000 0.0000 3.0000 0.0000" + LS
      + "4 0.0000 0.0000 0.0000 0.0000 2.0000 0.0000 0.0000" + LS
      ;

    assertEquals(exp, plot(ac));

    final String expRev = ""
      + "0 0.0000 2.0000 0.0000 0.0000 0.0000 0.0000 0.0000" + LS
      + "1 1.0000 0.0000 3.0000 0.0000 0.0000 0.0000 0.0000" + LS
      + "2 0.0000 0.0000 0.0000 0.0000 0.0000 3.0000 4.0000" + LS
      + "3 0.0000 0.0000 0.0000 0.0000 2.0000 0.0000 0.0000" + LS
      + "4 0.0000 0.0000 0.0000 1.0000 0.0000 0.0000 0.0000" + LS
      ;
    assertEquals(expRev, plot(ac.reverse()));

  }

  private String plot(final AllCounts ac) {
    final MemoryPrintStream ps = new MemoryPrintStream();
    ac.plot(ps.printStream());
    return ps.toString().replace("\t", " ");
  }
}
