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
public class MetagenomicsReportTemplateTest extends TestCase {
  public void testTemplate() throws IOException {
    final MetagenomicsReportTemplate template = new MetagenomicsReportTemplate();
    template.mSpeciesReport = "species report is here";
    template.mSpecies = "you asked for species to be run";
    template.mVersion = "A version file";
    template.mMapf = "mapf is all go";
    template.mMapfReport = "output from a hypthetical mapf run";
    template.mMapx = "more text";
    template.mMap = "checked";
    template.mMapxReport = "it all works";
    // this will error if there is there is a mismatch between template replace markers and class members
    template.fillTemplate();
  }
}
