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
public class NetworkStatsTest extends TestCase {

  public void test() throws IOException {
    NetworkStats ds = new NetworkStats("rubbish");
    StringBuffer sb = new StringBuffer();
    ds.addHeader(sb);
    assertEquals(0, sb.length());
    ds.addColumnLabelsTop(sb);
    assertEquals(0, sb.length());
    ds.addColumnLabelsBottom(sb);
    assertEquals(0, sb.length());
    ds.addColumnData(sb);
    assertEquals(0, sb.length());
  }

  public void test2() throws IOException {
    NetworkStats ds = new NetworkStats("eth0");
    StringBuffer sb = new StringBuffer();
    ds.addHeader(sb);
    assertEquals(0, sb.length());

    ds.addColumnLabelsTop(sb);
    if (sb.length() > 0) {
      assertTrue(sb.toString().contains("eth0"));

      sb.setLength(0);
      ds.addColumnData(sb);
      assertTrue(sb.length() > 0);
      assertTrue(sb.toString().contains("n/a"));
      sb.setLength(0);
      ds.addColumnData(sb);
      assertTrue(sb.length() > 0);
      assertFalse(sb.toString().contains("n/a"));
    }
  }

  public static Test suite() {
    return new TestSuite(NetworkStatsTest.class);
  }

  /**
   * Main to run from tests from command line.
   * @param args ignored.
   */
  public static void main(final String[] args) {
    junit.textui.TestRunner.run(NetworkStatsTest.class);
  }

}

