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

import java.util.HashMap;
import java.util.Map;

/**
 */
public class MetagenomicsReportTemplate extends ReportTemplate {
  String mVersion;
  String mSpecies;
  String mMap;
  String mSpeciesReport;
  String mMapx;
  String mMapxReport;
  String mMapf;
  String mMapfReport;

  @Override
  Map<String, String> makeReplacements() {
    final Map<String, String> replacements = new HashMap<>();
    replacements.put("VERSION", mVersion);
    replacements.put("SPECIES", mSpecies);
    replacements.put("MAP", mMap);
    replacements.put("SPECIES_REPORT", mSpeciesReport);
    replacements.put("MAPX", mMapx);
    replacements.put("MAPX_REPORT", mMapxReport);
    replacements.put("MAPF", mMapf);
    replacements.put("MAPF_REPORT", mMapfReport);
    return replacements;
  }
  MetagenomicsReportTemplate() {
    super("metagenomics_report.html");
  }
}
