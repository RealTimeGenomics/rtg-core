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
public class MonUtilsTest extends TestCase {

  public void test() throws IOException {
    double n = 0.13213;
    assertEquals("0.13", MonUtils.NF2.format(n));
    assertEquals("0.1", MonUtils.NF1.format(n));
    assertEquals("0", MonUtils.NF0.format(n));
    StringBuffer sb = new StringBuffer();
    MonUtils.pad(sb, "0", 3);
    assertEquals("  0", sb.toString());

    sb.setLength(0);
    MonUtils.padRight(sb, "0", 3);
    assertEquals("0  ", sb.toString());
  }


  public static Test suite() {
    return new TestSuite(MonUtilsTest.class);
  }

  /**
   * Main to run from tests from command line.
   * @param args ignored.
   */
  public static void main(final String[] args) {
    junit.textui.TestRunner.run(MonUtilsTest.class);
  }

}

