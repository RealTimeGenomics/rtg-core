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
