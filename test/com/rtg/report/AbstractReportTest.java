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

package com.rtg.report;

import java.io.File;

import com.rtg.util.diagnostic.Diagnostic;

import junit.framework.TestCase;

/**
 */
public abstract class AbstractReportTest extends TestCase {

  private String mHeadlessState = null;

  @Override
  public void setUp() throws Exception {
    super.setUp();
    Diagnostic.setLogStream();
    mHeadlessState = System.getProperty("java.awt.headless");
    System.setProperty("java.awt.headless", Boolean.TRUE.toString());
  }

  @Override
  public void tearDown() throws Exception {
    super.tearDown();
    if (mHeadlessState != null) {
      System.setProperty("java.awt.headless", mHeadlessState);
      mHeadlessState = null;
    } else {
      System.clearProperty("java.awt.headless");
    }
  }

  abstract Report getReport();

  public void testBadParams() throws Exception {
    final Report report = getReport();
    try {
      report.generateReport(ReportType.PDF, null, (File) null);
      fail("Accepted PDF report type.");
    } catch (UnsupportedOperationException uoe) {
      assertEquals("Only " + ReportType.HTML + " reports implemented so far", uoe.getMessage());
    }
    try {
      report.generateReport(ReportType.HTML, null, (File) null);
      fail("Accepted null directory.");
    } catch (NullPointerException npe) {
      assertNull(npe.getMessage());
    }
  }
}
