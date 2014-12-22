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
public class RecordStatsTest extends TestCase {

  private static final String LS = System.lineSeparator();

  public void test() throws IOException {
    StringBuffer sb = new StringBuffer();
    RecordStats ds = new RecordStats(sb, 10000);
    ds.addStats(new MonStats() {
        @Override
        public void addHeader(Appendable out) throws IOException {
          out.append("# Test").append(LS);
        }

        @Override
        public void addColumnLabelsTop(Appendable out) throws IOException {
          out.append("Top ");
        }

        @Override
        public void addColumnLabelsBottom(Appendable out) throws IOException {
          out.append("Bot ");
        }

        @Override
        public void addColumnData(Appendable out) throws IOException {
          out.append("Data");
        }

      });
    ds.addHeader();
    assertTrue(sb.toString(), sb.toString().startsWith("#"));
    assertEquals(2, sb.toString().split(LS).length);
    sb.setLength(0);
    ds.addColumnLabels();
    assertTrue(sb.toString(), sb.toString().startsWith("#"));
    assertEquals(3, sb.toString().split(LS).length);
    sb.setLength(0);
    ds.addColumnData();
    assertTrue(sb.toString(), sb.toString().contains("Data"));
    assertEquals(1, sb.toString().split(LS).length);
  }

  public static Test suite() {
    return new TestSuite(RecordStatsTest.class);
  }

  /**
   * Main to run from tests from command line.
   * @param args ignored.
   */
  public static void main(final String[] args) {
    junit.textui.TestRunner.run(RecordStatsTest.class);
  }

}

