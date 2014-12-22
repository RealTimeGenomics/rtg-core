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
public class DiskStatsTest extends TestCase {

  public void test() throws IOException {
    DiskStats ds = new DiskStats("rubbish");
    StringBuffer sb = new StringBuffer();
    ds.addHeader(sb);
    assertEquals(sb.toString(), 0, sb.length());
    sb.setLength(0);
    ds.addColumnData(sb);
    assertEquals(sb.toString(), 0, sb.length());
    sb.setLength(0);
    ds.addColumnLabelsTop(sb);
    assertEquals(sb.toString(), 0, sb.length());
    sb.setLength(0);
    ds.addColumnLabelsBottom(sb);
    assertEquals(sb.toString(), 0, sb.length());
  }

  public void test2() throws IOException {
    DiskStats ds = new DiskStats("sda");
    StringBuffer sb = new StringBuffer();
    ds.addHeader(sb);
    assertEquals(0, sb.length());

    ds.addColumnLabelsTop(sb);
    if (sb.length() > 0) {
      assertTrue(sb.toString().contains("sda"));

      sb.setLength(0);
      ds.addColumnData(sb);
      assertTrue(sb.length() > 0);
    }
  }

  public static Test suite() {
    return new TestSuite(DiskStatsTest.class);
  }

  /**
   * Main to run from tests from command line.
   * @param args ignored.
   */
  public static void main(final String[] args) {
    junit.textui.TestRunner.run(DiskStatsTest.class);
  }

}

