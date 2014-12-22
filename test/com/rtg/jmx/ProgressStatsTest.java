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

import com.rtg.util.cli.CommandLine;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

/**
 */
public class ProgressStatsTest extends TestCase {

  public void test() throws IOException {
    ProgressStats ds = new ProgressStats();
    StringBuffer sb = new StringBuffer();
    CommandLine.setCommandArgs("a", "simple", "command", "line");
    ds.addHeader(sb);
    assertTrue(sb.toString().startsWith("# Command line "));
    assertTrue(sb.toString().contains("a simple command line"));
  }

  public static Test suite() {
    return new TestSuite(ProgressStatsTest.class);
  }

  /**
   * Main to run from tests from command line.
   * @param args ignored.
   */
  public static void main(final String[] args) {
    junit.textui.TestRunner.run(ProgressStatsTest.class);
  }

}

