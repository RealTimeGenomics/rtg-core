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
package com.rtg.simulation;

/**
 * Test class
 */
public class OldReadNameParserTest extends AbstractSimulatedReadNameParserTest {

  @Override
  protected SimulatedReadNameParser getParser() {
    return new OldReadNameParser();
  }

  @Override
  protected String getMalformedReadName() {
    return "readname/adskfl://";
  }

  @Override
  protected String getWellFormedReadName(String tName, int tPos, boolean forward, int rLen) {
    final int newPos = forward ? tPos : tPos + rLen;
    return "read1:" + tName + ":" + newPos + (forward ? "" : "R") + ":S0:I0:D0";
  }

  public void testCigarScanning() {
    final OldReadNameParser p = new OldReadNameParser();
    assertFalse(p.setReadInfo("read9:tem:7R:S1:I2", 35));
    assertTrue(p.setReadInfo("read9:tem:7R:S1:I2:D3", 35));
    assertEquals(2, p.insertions());
    assertEquals(3, p.deletions());
    assertEquals(1, p.substitutions());
    assertEquals("read9", p.readName());
    assertEquals(9, p.readId());
    assertEquals(0, p.templateSet());
    assertFalse(p.forwardFrame());
    assertEquals(-28, p.templatePosition());
  }
}
