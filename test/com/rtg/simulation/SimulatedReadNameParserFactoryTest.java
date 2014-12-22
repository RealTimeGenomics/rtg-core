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

import junit.framework.TestCase;

/**
 * Test class
 */
public class SimulatedReadNameParserFactoryTest extends TestCase {

  public void testFactory() {
    SimulatedReadNameParser p = SimulatedReadNameParserFactory.getParser("read9 0/0/simulatedSequence1/9068/F/7.1X27.");
    assertNotNull(p);
    assertTrue(p.getClass().equals(NewReadNameParser.class));

    p = SimulatedReadNameParserFactory.getParser("read9/0/0/simulatedSequence1/9068/F/7.1X27.");
    assertNotNull(p);
    assertTrue(p.getClass().equals(NewReadNameParser.class));

    p = SimulatedReadNameParserFactory.getParser("cgread02:simulatedSequence1:30713R:S0:I0:D0");
    assertNotNull(p);
    assertTrue(p instanceof OldReadNameParser);

    p = SimulatedReadNameParserFactory.getParser("0 frag9/0/0/simulatedSequence1/9068/F/7.1X27.");
    assertNotNull(p);
    assertTrue(p instanceof NewestReadNameParser);

    p = SimulatedReadNameParserFactory.getParser("@CFTR.3.70s_282_148_0_1_0_0_4:0:0_4:0:0_10c/1");
    assertNotNull(p);
    assertTrue(p instanceof DwgsimReadNameParser);

    assertNull(SimulatedReadNameParserFactory.getParser("cgread2:si6"));
    assertNull(SimulatedReadNameParserFactory.getParser("cgread2/si6"));
    assertNull(SimulatedReadNameParserFactory.getParser("/////"));
    assertNull(SimulatedReadNameParserFactory.getParser("::::"));
    assertNull(SimulatedReadNameParserFactory.getParser(""));
  }
}
