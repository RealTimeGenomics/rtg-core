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

package com.rtg.variant.sv;

import static com.rtg.util.StringUtils.LS;

import com.rtg.util.TestUtils;
import com.rtg.util.io.MemoryPrintStream;

import junit.framework.TestCase;

/**
 */
public class AllCountsTest extends TestCase {

  public void test() {
    final AllCounts ac = new AllCounts(5);
    ac.globalIntegrity();
    assertEquals("AllCounts:5", ac.toString());
    TestUtils.containsAll(plot(ac), "0 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000", "4 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000");

    ac.properLeft().increment(0);
    ac.discordantLeft().increment(1, 2);
    ac.unmatedLeft().increment(2, 3);
    ac.properRight().increment(3);
    ac.discordantRight().increment(4, 2);
    ac.unmatedRight().increment(3, 3);
    ac.unpaired().increment(2, 4);
    final String exp = ""
      + "0 1.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000" + LS
      + "1 0.0000 2.0000 0.0000 0.0000 0.0000 0.0000 0.0000" + LS
      + "2 0.0000 0.0000 3.0000 0.0000 0.0000 0.0000 4.0000" + LS
      + "3 0.0000 0.0000 0.0000 1.0000 0.0000 3.0000 0.0000" + LS
      + "4 0.0000 0.0000 0.0000 0.0000 2.0000 0.0000 0.0000" + LS
      ;

    assertEquals(exp, plot(ac));

    final String expRev = ""
      + "0 0.0000 2.0000 0.0000 0.0000 0.0000 0.0000 0.0000" + LS
      + "1 1.0000 0.0000 3.0000 0.0000 0.0000 0.0000 0.0000" + LS
      + "2 0.0000 0.0000 0.0000 0.0000 0.0000 3.0000 4.0000" + LS
      + "3 0.0000 0.0000 0.0000 0.0000 2.0000 0.0000 0.0000" + LS
      + "4 0.0000 0.0000 0.0000 1.0000 0.0000 0.0000 0.0000" + LS
      ;
    assertEquals(expRev, plot(ac.reverse()));

  }

  private String plot(final AllCounts ac) {
    final MemoryPrintStream ps = new MemoryPrintStream();
    ac.plot(ps.printStream());
    return ps.toString().replace("\t", " ");
  }
}
