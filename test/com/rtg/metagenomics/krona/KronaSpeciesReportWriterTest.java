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
package com.rtg.metagenomics.krona;

import java.io.File;

import com.rtg.util.HtmlReportHelper;
import com.rtg.util.io.TestDirectory;

import junit.framework.TestCase;

/**
 */
public class KronaSpeciesReportWriterTest extends TestCase {

  public void testCopyResources() throws Exception {
    try (TestDirectory td = new TestDirectory()) {
      final HtmlReportHelper hrh = new HtmlReportHelper(td, "index");
      final KronaSpeciesReportWriter ksrw = new KronaSpeciesReportWriter(hrh);

      ksrw.copyResources();
      final File repDir = hrh.getResourcesDir();
      assertTrue(new File(repDir, "krona-2.0.js").exists());
      assertTrue(new File(repDir, "loading.gif").exists());
      assertTrue(new File(repDir, "hidden.png").exists());
    }
  }
}
