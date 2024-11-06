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
package com.rtg.position.output;

import static com.rtg.util.StringUtils.LS;

import java.io.IOException;
import java.io.StringWriter;

import com.rtg.mode.UnidirectionalFrame;

/**
 * Tests the corresponding class.
 *
 */
public class RawPositionWriterTest extends AbstractPositionWriterTest {

  public void test() throws IOException {
    final StringWriter out = new StringWriter();
    final StringWriter ambiguous = new StringWriter();
    final StringWriter unmapped = new StringWriter();
    final PositionWriter w = new RawPositionWriter(out, unmapped);
    assertEquals(0, w.score(), 1E-12);
    w.endQuery(1);
    w.write(getRegion(2, 1, 32, 32.0), 43, UnidirectionalFrame.FORWARD, 44, 0 /*unused*/);
    w.endQuery(43);
    w.write(getRegion(2, 1, 32, 32.0), 53, UnidirectionalFrame.FORWARD, 44, 0);
    w.write(getRegion(2, 1, 32, 32.0), 53, UnidirectionalFrame.FORWARD, 44, 0);
    w.endQuery(53);
    w.endQuery(63);
    assertEquals("1" + LS + "63" + LS, unmapped.toString());
    final String exp = ""
      + "#query-id\tquery-frame\tquery-start\tquery-end\tsubject-id\tsubject-frame\tsubject-start\tsubject-end" + LS
      + "43\t\t33\t64\t2\t\t2\t33" + LS
      + "53\t\t33\t64\t2\t\t2\t33" + LS
      + "53\t\t33\t64\t2\t\t2\t33" + LS
      ;
    assertEquals(exp, out.toString());
    assertEquals("", ambiguous.toString());
    assertEquals(96.0, w.score());
  }
}

