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
