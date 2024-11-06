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
public class NewReadNameParserTest extends AbstractSimulatedReadNameParserTest {


  @Override
  protected SimulatedReadNameParser getParser() {
    return new NewReadNameParser();
  }

  @Override
  protected String getMalformedReadName() {
    return "readname/adskfl://";
  }

  @Override
  protected String getWellFormedReadName(String tName, int tPos, boolean forward, int rLen) {
    return "read1 0/0/" + tName + "/" + tPos + "/" + (forward ? "F" : "R") + "/" + rLen + ".";
  }

  public void testCigarScanning() {
    final NewReadNameParser p = new NewReadNameParser();

    assertFalse(p.setReadInfo("read9 0/0/simulatedSequence1/9068/F", 35));
    assertTrue(p.setReadInfo("read9 0/0/simulatedSequence1/9068/F/7.1X27./Right", 35));
    assertEquals(0, p.insertions());
    assertEquals(0, p.deletions());
    assertEquals(1, p.substitutions());
    assertEquals("7=1X27=", p.cigar());

    p.setReadInfo("read1 0/8/simulatedSequence9/65092/R/13.2X10.1X5.1X2.1X/Right", 35);
    assertEquals(0, p.insertions());
    assertEquals(0, p.deletions());
    assertEquals(5, p.substitutions());
    assertEquals("13=2X10=1X5=1X2=1X", p.cigar());

    p.setReadInfo("read73/0/0/simulatedSequence1/106/F/14.1X15.", 30);
    assertEquals(106, p.templatePosition());
    assertTrue(p.forwardFrame());
    assertEquals("14=1X15=", p.cigar());

    assertTrue(p.setReadInfo("read1 0/8/simulatedSequence9/65092/R/0X13.20D1I12.", 35));
    assertEquals(1, p.insertions());
    assertEquals(20, p.deletions());
    assertEquals(0, p.substitutions());
    assertEquals("0X13=20D1I12=", p.cigar());
    assertEquals("read1", p.readName());
    assertEquals(1, p.readId());
    assertEquals(0, p.templateSet());
  }

}
