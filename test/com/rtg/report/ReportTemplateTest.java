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
