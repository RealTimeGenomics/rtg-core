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
package com.rtg.util;

import java.io.File;

import com.rtg.util.io.TestDirectory;

import junit.framework.TestCase;

/**
 */
public class HtmlReportHelperTest extends TestCase {

  public void testBits() throws Exception {
    try (TestDirectory dir = new TestDirectory()) {
      final HtmlReportHelper hrh = new HtmlReportHelper(dir, "report");
      assertEquals("report_files", hrh.getResourcesDirName());
      assertEquals("report.html", hrh.getReportFile().getName());
      assertEquals(dir, hrh.getReportFile().getParentFile());
      assertEquals(dir, hrh.getBaseDir());

      hrh.copyResources("com/rtg/util/resources/krona.xml");
      final File repDir = new File(dir, hrh.getResourcesDirName());
      assertEquals(hrh.getResourcesDir(), repDir);
      assertTrue(repDir.exists());
      final File f = new File(repDir, "krona.xml");
      assertTrue(f.exists());
    }
  }
}
