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
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import com.rtg.util.Resources;
import com.rtg.util.StringUtils;

/**
 * Manages the replacement of markers within a template file with the fields marked with the
 * Searches for tags within the template that look like <code>__SOME_TAG_NAME__</code> and replaces
 * them with a value as specified by <code>makeReplacements</code>
 * Note: there is a find bugs exclusion for Unread fields that would otherwise be triggered by these
 */
abstract class ReportTemplate {
  static final Pattern REPLACE = Pattern.compile("__([A-Z]|_[A-Z])+__");
  final String mTemplateFile;
  ReportTemplate(String templateFile) {
    mTemplateFile = templateFile;
  }

  /**
   * Find the fields the members that should be used and construct a map of tag name to value
   * @return a map from tag name (without leading/trailing underscores) to replacement value
   */
  abstract Map<String, String> makeReplacements();

  /**
   * replaces tags in the template with the <code>TemplateReplacement</code> values
   * @return the populated template as a string
   * @throws IOException if template reading fails
   */
  public String fillTemplate() throws IOException {
    final Map<String, String> replacements = makeReplacements();
    final Map<String, Boolean> wasPopulated = new HashMap<>();
    for (Map.Entry<String, String> entry : replacements.entrySet()) {
      if (entry.getValue() == null) {
        throw new IllegalArgumentException(entry.getKey() + " is not set");
      }
    }

    for (String key : replacements.keySet()) {
      wasPopulated.put(key, false);
    }
    final StringBuilder body = new StringBuilder();
    final String templateFile = ReportUtils.TEMPLATE_DIR + "/" + mTemplateFile;
    try (BufferedReader template = new BufferedReader(new InputStreamReader(Resources.getResourceAsStream(SpeciesReport.class, templateFile)))
    ) {
      String line;
      while ((line = template.readLine()) != null) {
        final Matcher matcher = REPLACE.matcher(line);
        int start = 0;
        while (matcher.find(start)) {
          final String match = matcher.group();
          final String key = match.substring(2, match.length() - 2);
          if (!replacements.keySet().contains(key)) {
            throw new IllegalArgumentException("Template contained unsupported replacement tag: " + key);
          }
          start = matcher.end();
        }
        for (Map.Entry<String, String> entry : replacements.entrySet()) {
          final String replaceKey = "__" + entry.getKey() + "__";
          if (line.contains(replaceKey)) {
            wasPopulated.put(entry.getKey(), true);
            line = line.replaceAll(replaceKey, entry.getValue().replaceAll("\\\\", "\\\\\\\\"));
          }
        }
        body.append(line).append(StringUtils.LS);
      }
    }
    for (Map.Entry<String, Boolean> entry : wasPopulated.entrySet()) {
      if (!entry.getValue()) {
        throw new IllegalArgumentException("Template: " + getClass().getName() + " had nowhere to put replacement " + entry.getKey());
      }
    }
    return body.toString();
  }
}
