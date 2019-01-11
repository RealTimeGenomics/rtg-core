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

import java.io.IOException;

import junit.framework.TestCase;

/**
 */
public class MapReportTemplateTest extends TestCase {
  // Basically just testing that all the parameters are present in the template
  public void testTemplate() throws IOException {
    final MapReportTemplate template = new MapReportTemplate();
    template.mFullReportTitle = "fred";
    template.mSummaryFile = "summary";
    template.mSummaryFileTitle = "super original title";
    template.mParameterSummary = " I am some parameters";
    template.mFullReport = "full report lives here";
    template.fillTemplate();
  }
}
