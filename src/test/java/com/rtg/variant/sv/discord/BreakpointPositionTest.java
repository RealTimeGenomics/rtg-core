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

package com.rtg.variant.sv.discord;

import junit.framework.TestCase;

/**
 */
public class BreakpointPositionTest extends TestCase {

  public void test() {
    final BreakpointPosition pos = new BreakpointPosition(1, 2, 3, 4);
    assertEquals(1, pos.lo());
    assertEquals(2, pos.position());
    assertEquals(3, pos.hi());
    assertEquals(4, pos.positionAlt());
    assertEquals("BreakpointPosition [lo=1, position=2, hi=3, positionY=4]", pos.toString());
  }
}
