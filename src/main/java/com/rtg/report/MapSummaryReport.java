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
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Locale;

import com.reeltwo.jumble.annotations.TestClass;
import com.reeltwo.plot.ui.ImageWriter;
import com.rtg.ngs.NgsParams;
import com.rtg.util.HtmlReportHelper;
import com.rtg.util.IntegerOrPercentage;
import com.rtg.util.StringUtils;
import com.rtg.util.diagnostic.Diagnostic;

/**
 * Write params used in mapping run to a report
 */
@TestClass("com.rtg.report.MapReportTest")
public class MapSummaryReport implements Report {

  static {
    System.setProperty("java.awt.headless", "true");
  }

  static final String REPORT_NAME = "index";

  /** Comment contents for summary delimiter */
  public static final String SUMMARY_STRING = "PARAMETER SUMMARY";
  /** Used to delimit the start of the summary block for use in combined report */
  public static final String START_PARAMS_SUMMARY = "<!--" + SUMMARY_STRING + "-->";
  /** Comment contents for summary end delimiter */
  public static final String END_SUMMARY_STRING = "END PARAMETER SUMMARY";
  /** End of summary block delimiter */
  public static final String END_PARAMS_SUMMARY = "<!--" + END_SUMMARY_STRING + "-->";
  private List<File> mSampleFiles;
  NgsParams mParams;
  private String mTitle = "Mapping Report";
  private List<String> mCommandLines = Collections.emptyList();

  /**
   * Constructor for MapReport
   */
  public MapSummaryReport() {
  }

  public void setTitle(String title) {
    mTitle = title;
  }

  public void setSampleFiles(List<File> sampleFiles) {
    mSampleFiles = sampleFiles;
  }

  public void setCommandLines(List<String> commandLines) {
    mCommandLines = commandLines;
  }

  public void setParams(NgsParams params) {
    mParams = params;
  }

  @Override
  public void generateReport(ReportType type, File outputDir, File... inputDirs) throws IOException {
    if (type != ReportType.HTML) {
      throw new UnsupportedOperationException("Only " + ReportType.HTML + " reports implemented so far");
    }
    if (inputDirs == null || outputDir == null) {
      throw new NullPointerException();
    }
    if (inputDirs.length == 0) {
      throw new IllegalArgumentException("No input directories given.");
    }
    for (File dir : inputDirs) {
      if (dir == null) {
        throw new NullPointerException();
      }
    }
    final String m = ImageWriter.isImageWritingEnabled();
    if (m != null) {
      Diagnostic.warning("Skipping HTML report output, host OS is not correctly configured: " + m);
    } else {
      writeReport(outputDir, inputDirs);
    }
  }

  private void writeReport(File outputDir, File... inputDirs) throws IOException {
    final HtmlReportHelper helper = new HtmlReportHelper(outputDir, REPORT_NAME);
    helper.copyResources(ReportUtils.TEMPLATE_DIR + "/rtg.css"); //copy resources up here to create the report files sub dir as well
    final StringBuilder body = new StringBuilder();

    reportContent(helper, body, inputDirs);

    // write html
    final HashMap<String, String> replacements = new HashMap<>();
    replacements.put("__TITLE__", mTitle);
    replacements.put("__BODY_TEXT__", body.toString());
    replacements.put("__RESOURCE_DIR__", helper.getResourcesDirName() + "/");
    ReportUtils.writeHtml(ReportUtils.TEMPLATE_DIR + "/default.html", helper.getReportFile(), replacements);
  }
  List<String> getCommandLines() {
    return mCommandLines;
  }

  void reportContent(HtmlReportHelper helper, StringBuilder body, File[] inputDirs) throws IOException {
    body.append(StringUtils.LS).append(START_PARAMS_SUMMARY).append(StringUtils.LS);
    writeParamsBlock(body);
    final List<String> commandLines = getCommandLines();
    if (commandLines.size() == 1) {
        body.append("<p><strong>Command line: ").append("</strong>").append(commandLines.get(0)).append("</p>");
      } else {
      body.append("<p><strong>Command lines").append("</strong></p>");
      body.append("<ul>");
      for (String cl : commandLines) {
        body.append("<li>");
        body.append(cl)
          .append(StringUtils.LS);
        body.append("</li>");
      }
      body.append("</ul>");
    }
    body.append(StringUtils.LS).append(END_PARAMS_SUMMARY).append(StringUtils.LS);

  }

  private void writeParamsBlock(StringBuilder body) throws IOException {
    if (mSampleFiles != null) {
      body.append("<strong>Sample:</strong> ");
      String join = "";
      for (File f : mSampleFiles) {
        if (f != null) {
          body.append(join);
          body.append(f.getPath());
          join = ", ";
        } else {
          throw new IllegalArgumentException();
        }
      }
      body.append("<br/>").append(StringUtils.LS);
    }
    if (mParams != null) {
      final String readMe =  mParams.searchParams().reader().getReadMe();
      final File template = mParams.searchParams().reader().path();
      if (readMe != null) {
        body.append("<strong>Database:</strong> ");
        body.append(readMe);
        body.append("<br/>");
        body.append(StringUtils.LS);
      } else if (template != null) {
        body.append("<strong>Database:</strong> ");
        body.append(template.getPath());
        body.append("<br/>");
        body.append(StringUtils.LS);
      }
      final File output = mParams.directory();
      if (output != null) {
        body.append("<strong>Results directory:</strong> ")
            .append(output.getPath())
            .append("<br/>")
            .append(StringUtils.LS);
      }

      sensitivityReport(body);
      reportingReport(body);
    }
  }

  void reportingReport(StringBuilder body) {
    final int topN = mParams.outputParams().maxTopResults();
    final String reportingFormatString = ""
                                         + "<p><strong>Reporting (SAM/BAM records):</strong> Top %d positions, uniquely mapping reads NH = 1, ambiguously"
                                         + " mapping reads up to %d positions, NH >1. Reads with more than %d hits are reported as unmapped (XC=C)."
                                         + " Reads with mappings below alignment threshold reported as XC=D. </p>";
    final String reportingString = String.format(Locale.ROOT, reportingFormatString, topN, topN, topN);
    body.append(reportingString);
  }

  void sensitivityReport(StringBuilder body) {
    int wordSize;
    try {
      wordSize = mParams.maskParams().getWordSize();
    } catch (UnsupportedOperationException e) {
      wordSize = -1;
    }
    boolean started = false;
    body.append("<p><strong>Sensitivity settings:</strong> ");
    if (wordSize >= 0) {
      body.append("Word size=")
          .append(wordSize);
      started = true;
    }
    final int stepSize = mParams.stepSize();
    if (stepSize > 0) {
      if (started) {
        body.append(", step size=");
      } else {
        body.append("Step size=");
        started = true;
      }
      body.append(stepSize);
    }
    if (started) {
        body.append(", alignment threshold=");
    } else {
      body.append("Alignment threshold=");
    }
    final IntegerOrPercentage matedAs = mParams.outputParams().filter().unmatedMaxMismatches();
    final IntegerOrPercentage unmatedAs = mParams.outputParams().filter().unmatedMaxMismatches();
    if (mParams.paired() && !unmatedAs.equals(matedAs)) {
      body.append(matedAs).append(" mated/").append(unmatedAs).append(" unmated");
    } else {

      body.append(unmatedAs);
      if (unmatedAs.isPercentage()) {
        body.append(" (").append(unmatedAs).append(" mismatch allowed over length of read)");
      }
    }
    body.append(".</p>");
  }

}
