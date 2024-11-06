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
package com.rtg.alignment;

import com.rtg.mode.DnaUtils;
import com.rtg.util.MaxShiftUtils;
import com.rtg.util.TestUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.MemoryPrintStream;

import junit.framework.TestCase;

/**
 */
public class LoggingOnlyEditDistanceTest extends TestCase {

  /**
   * To test alignment speed, set repeatCount above to 100,000 then run this
   * test.
   */
  public void testEmpty() {

    MemoryPrintStream ms = new MemoryPrintStream();
    Diagnostic.setLogStream(ms.printStream());
    final byte[] read = new byte[0];
    final byte[] template = new byte[0];
    UnidirectionalEditDistance lo = new LoggingOnlyEditDistance();
    int[] res = lo.calculateEditDistance(read, 0, template, 0, 0, MaxShiftUtils.calculateDefaultMaxShift(read.length), true);
    assertNull(res);
    // Logging only Test does just that?
    TestUtils.containsAll(ms.toString(), "LoggingOnlyEditDistance  1:", "crc=0", "read.length=0", "rlen=0", "zeroStart=0", "template.length=0",
        "maxScore=0", "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN");

    lo = new LoggingOnlyEditDistance(3, "yo");
    ms = new MemoryPrintStream();
    Diagnostic.setLogStream(ms.printStream());
    res = lo.calculateEditDistance(DnaUtils.encodeString("AAAGCCCG"), 7, DnaUtils.encodeString("TTAAAGGCCGTT"), 1, 2, MaxShiftUtils.calculateDefaultMaxShift(read.length), true);
    assertNull(res);
    TestUtils.containsAll(ms.toString(), "LoggingOnlyEditDistance yo 1:", "crc=3036913", "read.length=8", "rlen=7", "zeroStart=1", "template.length=12",
        "maxScore=2", "ED: AAAGCCCG TTAAAGGCCGTTNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN 20");

    ms = new MemoryPrintStream();
    Diagnostic.setLogStream(ms.printStream());
    res = lo.calculateEditDistance(DnaUtils.encodeString("AAAGCCCG"), 7, DnaUtils.encodeString("TTTTTGGGGGAAAAACCCCCTTAAAGGCCGTT"), 21, 2, MaxShiftUtils.calculateDefaultMaxShift(read.length), true);
    assertNull(res);
    TestUtils.containsAll(ms.toString(), "ED: AAAGCCCG TTTTGGGGGAAAAACCCCCTTAAAGGCCGTTNNNNNNNNNNNNNNNN 20");
    Diagnostic.setLogStream();
  }



}
