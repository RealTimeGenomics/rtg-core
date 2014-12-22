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

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.util.Map;

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.util.Resources;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.LineWriter;

/**
 * Helper utils for writing reports
 */
@TestClass("com.rtg.report.MapReportTest")
public final class ReportUtils {
  /** location of resources */
  public static final String TEMPLATE_DIR = "com/rtg/report/resources";
  private ReportUtils() { }


  /**
   * Truncates string at given length, appends <code>"..."</code> if string was shortened
   * @param string string to truncate
   * @param length length to truncate at
   * @return if string length is less than or equal to given length then the given string, otherwise the given string truncated at given length with {@code ...} appended
   */
  public static String truncate(String string, int length) {
    return string.length() > length ? string.substring(0, length) + "..." : string;
  }

  static String[] resourceArray(String... resources) {
    final String[] result = new String[resources.length];
    for (int i = 0; i < resources.length; i++) {
      final String r = resources[i];
      result[i] = TEMPLATE_DIR + "/" + r;
    }
    return result;
  }

  /**
   * writes given template resource to given output file replacing things in <code>replacements</code> when found
   * @param templateFile name ot resource containing template
   * @param outputFile output file
   * @param replacements mapping between placeholder and final text
   * @throws IOException if an IO error occurs
   */
  public static void writeHtml(String templateFile, File outputFile, Map<String, String> replacements) throws IOException {
    try (final OutputStream os = FileUtils.createOutputStream(outputFile)) {
      writeHtml(templateFile, os, replacements);
    }
  }

  /**
   * writes given template resource to given output stream replacing things in <code>replacements</code> when found
   * @param templateFile name ot resource containing template
   * @param output output stream
   * @param replacements mapping between placeholder and final text
   * @throws IOException if an IO error occurs
   */
  public static void writeHtml(String templateFile, OutputStream output, Map<String, String> replacements) throws IOException {
    final InputStream is = Resources.getResourceAsStream(templateFile);

    try (BufferedReader br = new BufferedReader(new InputStreamReader(is))) {
      final LineWriter pw = new LineWriter(new OutputStreamWriter(output));
      fillTemplate(replacements, br, pw);
    }
  }

  static void fillTemplate(Map<String, String> replacements, BufferedReader br, LineWriter pw) throws IOException {
    try {
      String line;
      while ((line = br.readLine()) != null) {
        for (Map.Entry<String, String> entry : replacements.entrySet()) {
          line = line.replace(entry.getKey(), entry.getValue());
        }
        pw.writeln(line);
      }
    } finally {
      pw.close();
    }
  }
}
