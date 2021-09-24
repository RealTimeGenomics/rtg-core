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
class MapReportTemplate extends ReportTemplate {
  String mParameterSummary;
  String mSummaryFile;
  String mSummaryFileTitle;
  String mFullReport;
  String mFullReportTitle;

  @Override
  Map<String, String> makeReplacements() {
    final Map<String, String> replacements = new HashMap<>();
    replacements.put("PARAMETER_SUMMARY", mParameterSummary);
    replacements.put("SUMMARY_FILE", mSummaryFile);
    replacements.put("SUMMARY_FILE_TITLE", mSummaryFileTitle);
    replacements.put("FULL_REPORT", mFullReport);
    replacements.put("FULL_NAME", mFullReportTitle);
    return replacements;
  }

  MapReportTemplate() {
    super("map_report.html");
  }
}
