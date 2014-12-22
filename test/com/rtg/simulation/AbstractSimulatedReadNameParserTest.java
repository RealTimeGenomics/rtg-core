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
 * Test class.
 */
public abstract class AbstractSimulatedReadNameParserTest extends TestCase {

  protected abstract SimulatedReadNameParser getParser();

  protected abstract String getWellFormedReadName(String tName, int tPos, boolean forward, int rLen);

  protected abstract String getMalformedReadName();

  public void test() throws Exception {
    SimulatedReadNameParser p = getParser();
    assertFalse(p.setReadInfo(getMalformedReadName(), 10));

    String rname = getWellFormedReadName("mytemplate", 100, true, 10);
    assertTrue(p.setReadInfo(rname, 10));
    assertEquals("mytemplate", p.templateName());
    assertEquals(100, p.templatePosition());
    assertTrue(p.forwardFrame());
  }
}
