/*
 * Copyright (c) 2018. Real Time Genomics Limited.
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

package com.rtg.report;

import java.io.File;
import java.io.IOException;

import com.rtg.AbstractTest;

/**
 */
public abstract class AbstractReportTest extends AbstractTest {

  private String mHeadlessState = null;

  @Override
  public void setUp() throws IOException {
    super.setUp();
    mHeadlessState = System.getProperty("java.awt.headless");
    System.setProperty("java.awt.headless", Boolean.TRUE.toString());
  }

  @Override
  public void tearDown() throws IOException {
    if (mHeadlessState != null) {
      System.setProperty("java.awt.headless", mHeadlessState);
      mHeadlessState = null;
    } else {
      System.clearProperty("java.awt.headless");
    }
    super.tearDown();
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
