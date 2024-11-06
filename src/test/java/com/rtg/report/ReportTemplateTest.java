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

import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import java.util.regex.Matcher;

import com.rtg.util.StringUtils;

import junit.framework.TestCase;

/**
 */
public class ReportTemplateTest extends TestCase  {
  private static class SillyReportTemplate extends ReportTemplate {
    SillyReportTemplate() {
      super("test.html");
    }
    String mReplace;
    String mSecondReplacement;
    @Override
    Map<String, String> makeReplacements() {
      final Map<String, String> replacements = new HashMap<>();
      replacements.put("REPLACE", mReplace);
      replacements.put("SECOND_REPLACEMENT", mSecondReplacement);
      return replacements;
    }
  }
  public void testReportTemplateSuccess() throws IOException {
    final SillyReportTemplate t = new SillyReportTemplate();
    t.mReplace = "replacement";
    t.mSecondReplacement = "second";
    final String s = t.fillTemplate();
    assertEquals("test" + StringUtils.LS + "replacement before second after" + StringUtils.LS, s);
  }

  public void testUnset() throws IOException {
    final SillyReportTemplate t = new SillyReportTemplate();
    t.mReplace = "replacement";
    try {
    t.fillTemplate();
      fail("Expected an Illegal argument exception because mSecond wasn't set");
    } catch (IllegalArgumentException e) {
      // Expected
    }
  }

  private static class NoSecondReportTemplate extends ReportTemplate {
    NoSecondReportTemplate() {
      super("test.html");
    }
    String mReplace;
    @Override
    Map<String, String> makeReplacements() {
      final Map<String, String> replacements = new HashMap<>();
      replacements.put("REPLACE", mReplace);
      return replacements;
    }

  }
  public void testMissingField() throws IOException {
    final NoSecondReportTemplate t = new NoSecondReportTemplate();
    t.mReplace = "replacement";
    try {
      t.fillTemplate();
      fail("Expected an IllegalArgumentException because NoSecond report doesn't define t.mSecondReplacement");
    } catch (IllegalArgumentException e) {
      assertTrue(e.getMessage(), e.getMessage().contains("Template contained unsupported replacement tag: SECOND_REPLACEMENT"));
    }
  }
  private static class NoReplacementReportTemplate extends ReportTemplate {
    NoReplacementReportTemplate() {
      super("test.html");
    }
    String mSecondReplacement;
    @Override
    Map<String, String> makeReplacements() {
      final Map<String, String> replacements = new HashMap<>();
      replacements.put("SECOND_REPLACEMENT", mSecondReplacement);
      return replacements;
    }

  }
  public void testMissingField2() throws IOException {
    final NoReplacementReportTemplate t = new NoReplacementReportTemplate();
    t.mSecondReplacement = "replacement";
    try {
      t.fillTemplate();
      fail("Expected an IllegalArgumentException because NoSecond report doesn't define t.mReplace");
    } catch (IllegalArgumentException e) {
      assertTrue(e.getMessage(), e.getMessage().contains("Template contained unsupported replacement tag: REPLACE"));
    }
  }

  private static class BonusReportTemplate extends ReportTemplate {
    BonusReportTemplate() {
      super("test.html");
    }
    String mReplace;
    String mSecondReplacement;
    String mBonus;
    @Override
    Map<String, String> makeReplacements() {
      final Map<String, String> replacements = new HashMap<>();
      replacements.put("REPLACE", mReplace);
      replacements.put("SECOND_REPLACEMENT", mSecondReplacement);
      replacements.put("BONUS", mBonus);
      return replacements;
    }

  }

  public void testBonusField() throws IOException {
    final BonusReportTemplate t = new BonusReportTemplate();
    t.mReplace = "replacement";
    t.mSecondReplacement = "second replacement that has spaces";
    t.mBonus = "bonus";
    try {
      t.fillTemplate();
      fail("Expected an IllegalArgumentException because __BONUS__ doesn't appear in index.html");
    } catch (IllegalArgumentException e) {
      assertTrue(e.getMessage(), e.getMessage().contains("had nowhere to put replacement BONUS"));
    }
  }

  public void testPattern() {
    final Matcher matcher = ReportTemplate.REPLACE.matcher("__ASDF__");
    matcher.find();
    assertEquals("__ASDF__", matcher.group());
    final Matcher matcher2 = ReportTemplate.REPLACE.matcher(" some arbitrary text __ASDF_GASD__ that fills __SECOND_ONE__ an entire line");
    matcher2.find();
    assertEquals("__ASDF_GASD__", matcher2.group());
    matcher2.find(matcher2.end());
    assertEquals("__SECOND_ONE__", matcher2.group());
  }
}
