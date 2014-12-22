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
package com.rtg.util.diagnostic;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

/**
 * Tests the corresponding class.
 *
 */
public class ParallelProgressTest extends TestCase {

  /**
   */
  public ParallelProgressTest(final String name) {
    super(name);
  }
  public static Test suite() {
    return new TestSuite(ParallelProgressTest.class);
  }
  /**
   * Main to run from tests from command line.
   * @param args ignored.
   */
  public static void main(final String[] args) {
    junit.textui.TestRunner.run(suite());
  }

  @Override
  public void tearDown() {
    Diagnostic.setLogStream();
  }

  public void test() throws Exception {
    final ByteArrayOutputStream bos = new ByteArrayOutputStream();
    try {
      try (PrintStream ps = new PrintStream(bos)) {
        Diagnostic.setProgressStream(ps);
        final ParallelProgress pp = new ParallelProgress("test");
        pp.updateProgress(0);
        pp.updateProgress(1);
        pp.close();
      }
    } finally {
      bos.close();
    }
    final String t = bos.toString();
    //System.err.println(t);
    assertTrue(t.contains("Starting: test"));
    assertTrue(t.contains("Processed 0% of test"));
    assertTrue(t.contains("Processed 1% of test"));
    assertTrue(t.contains("Finished: test"));
  }

 }

