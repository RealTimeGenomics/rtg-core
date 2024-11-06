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
    for (Map.Entry<String, String> entry : replacements.entrySet()) {
      if (entry.getValue() == null) {
        throw new IllegalArgumentException(entry.getKey() + " is not set");
      }
    }

    final Map<String, Boolean> wasPopulated = new HashMap<>(replacements.size());
    for (String key : replacements.keySet()) {
      wasPopulated.put(key, Boolean.FALSE);
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
          if (!replacements.containsKey(key)) {
            throw new IllegalArgumentException("Template contained unsupported replacement tag: " + key);
          }
          start = matcher.end();
        }
        for (Map.Entry<String, String> entry : replacements.entrySet()) {
          final String replaceKey = "__" + entry.getKey() + "__";
          if (line.contains(replaceKey)) {
            wasPopulated.put(entry.getKey(), Boolean.TRUE);
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
