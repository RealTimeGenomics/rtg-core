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
package com.rtg.jmx;

import java.io.IOException;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

/**
 */
public class MBeanStatsTest extends TestCase {

  private static final String LS = System.lineSeparator();

  public void test() throws IOException {
    MBeanStats ds = new MBeanStats();
    StringBuffer sb = new StringBuffer();
    ds.addHeader(sb);
    assertEquals(8, sb.toString().split(LS).length);
  }

  public static Test suite() {
    return new TestSuite(MBeanStatsTest.class);
  }

  /**
   * Main to run from tests from command line.
   * @param args ignored.
   */
  public static void main(final String[] args) {
    junit.textui.TestRunner.run(MBeanStatsTest.class);
  }

}

