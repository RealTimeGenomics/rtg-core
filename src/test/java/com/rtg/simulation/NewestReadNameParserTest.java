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
public class NewestReadNameParserTest extends AbstractSimulatedReadNameParserTest {


  @Override
  protected SimulatedReadNameParser getParser() {
    return new NewestReadNameParser();
  }

  @Override
  protected String getMalformedReadName() {
    return "readname/adskfl://";
  }

  @Override
  protected String getWellFormedReadName(String tName, int tPos, boolean forward, int rLen) {
    return "12 frag1/0/0/" + tName + "/" + tPos + "/" + (forward ? "F" : "R") + "/" + rLen + ".";
  }

  public void testCigarScanning() {
    final NewestReadNameParser p = new NewestReadNameParser();

    assertFalse(p.setReadInfo("12 frag9/0/0/simulatedSequence1/9068/F", 35));
    assertFalse(p.setReadInfo("read9 0/0/simulatedSequence1/9068/F", 35));
    assertTrue(p.setReadInfo("12 frag9/0/0/simulatedSequence1/9068/F/7.1X27./Right", 35));
    assertEquals(0, p.insertions());
    assertEquals(0, p.deletions());
    assertEquals(1, p.substitutions());
    assertEquals("7=1X27=", p.cigar());
    assertEquals("frag9", p.readName());
    assertEquals(12, p.readId());
  }

}
