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

package com.rtg.variant.sv.discord.pattern;

import junit.framework.TestCase;

/**
 *         Date: 15/03/12
 *         Time: 3:29 PM
 */
public class BreakpointAltTest extends TestCase {
  void check(String orig, String remoteChr, int remotePos, boolean localUp, boolean remoteUp) {
    BreakpointAlt b = new BreakpointAlt(orig);
    assertEquals(remoteChr, b.getRemoteChr());
    assertEquals(remotePos, b.getRemotePos());
    assertEquals(localUp, b.isLocalUp());
    assertEquals(remoteUp, b.isRemoteUp());
  }
  public void testBreakpointAlt() {
    check("A[foo:221[", "foo", 221, true, false);
    check("A]foo:221]", "foo", 221, true, true);
    check("[foo:221[A", "foo", 221, false, false);
    check("]foo:221]A", "foo", 221, false, true);
    check("]foo:261]A", "foo", 261, false, true);
    check("]bar:261]A", "bar", 261, false, true);

  }
}
