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
package com.rtg;


import java.io.ByteArrayOutputStream;
import java.util.Locale;
import java.util.ResourceBundle;

import com.rtg.util.StringUtils;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

/**
 * Test for the <code>License</code> class
 */
public class LicenseCommandTest extends TestCase {

  private static final ResourceBundle USAGE_RESOURCE = ResourceBundle.getBundle("com.rtg.Usage", Locale.getDefault());

  public void testModuleList() {
    final ByteArrayOutputStream baos = new ByteArrayOutputStream();
    assertEquals(0, LicenseCommand.mainInit(baos, CoreCommand.INFO));

    final String result = baos.toString();
    assertTrue(result.contains("License: "));
    assertTrue(result.contains("Licensed to: "));
    int longestUsageLength = "Command name".length();
    // Get longest string lengths for use below in pretty-printing.
    for (Command module : CoreCommand.INFO.commands()) {
      if (module.getCommandName().length() > longestUsageLength && LicenseCommand.showModule(module)) {
        longestUsageLength = module.getCommandName().length();
      }
    }
    final StringBuilder extraSpaces = new StringBuilder();
    for (int i = 0; i < longestUsageLength - "Command name".length(); i++) {
      extraSpaces.append(" ");
    }
    assertTrue(result.contains("\tCommand name" + extraSpaces + " \tLicensed?         Release Level" + StringUtils.LS));
    for (Command module : CoreCommand.INFO.commands()) {
      if (module.isLicensed() || !module.isHidden()) {
        assertTrue(module.getCommandName() + " not found", contains(result, module.getCommandName(), module.isLicensed(), module.getReleaseLevel(), longestUsageLength));
        assertTrue(module.getCategory().toString() + " not found", result.contains(StringUtils.LS + StringUtils.LS + USAGE_RESOURCE.getString(module.getCategory().toString() + "_CAT_DESC") + ":" + StringUtils.LS));
      } else {
        assertFalse(module.getCommandName() + " found!", contains(result, module.getCommandName(), false, module.getReleaseLevel(), longestUsageLength));
      }
    }
  }

  private static boolean contains(String haystack, String needle, boolean licensed, ReleaseLevel releaseLevel, int longestUsageLength) {
    final StringBuilder extraSpaces = new StringBuilder();
    for (int i = 0; i < longestUsageLength - needle.length(); i++) {
      extraSpaces.append(" ");
    }
    final String rLevel;
    if (releaseLevel == ReleaseLevel.GA) {
      rLevel = releaseLevel.toString();
    } else {
      rLevel = releaseLevel.toString().toLowerCase(Locale.getDefault());
    }
    return haystack.contains("\t" + needle.toLowerCase(Locale.getDefault()) + extraSpaces + " \t" + LicenseCommand.padTo(licensed ? "Licensed" : "Unlicensed", 17) + " " + rLevel);
  }

  public void testShowModule() {
    for (Command mod : CoreCommand.INFO.commands()) {
      switch (mod.getReleaseLevel()) {
        case ALPHA:
        case BETA:
          if (mod.isLicensed()) {
            assertTrue(LicenseCommand.showModule(mod));
          } else {
            assertFalse(LicenseCommand.showModule(mod));
          }
          break;
        case GA:
          assertTrue(LicenseCommand.showModule(mod));
          break;
        default:
          break;
      }
    }
  }


  public void testLongestLengthModule() {
    assertEquals(3, LicenseCommand.getLongestLengthModule(new Command[0], "abc"));
    assertEquals(6, LicenseCommand.getLongestLengthModule(new Command[]{CoreCommand.FORMAT.module()}, "abc"));
  }

  public static Test suite() {
    return new TestSuite(LicenseCommandTest.class);
  }

  public static void main(final String[] args) {
    junit.textui.TestRunner.run(new TestSuite(LicenseCommandTest.class));
  }
}
