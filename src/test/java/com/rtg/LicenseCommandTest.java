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
package com.rtg;


import java.io.ByteArrayOutputStream;
import java.util.Locale;

import com.rtg.launcher.globals.GlobalFlags;
import com.rtg.util.StringUtils;

import junit.framework.TestCase;

/**
 * Test for the <code>License</code> class
 */
public class LicenseCommandTest extends TestCase {

  public void testModuleList() {
    GlobalFlags.resetAccessedStatus();
    final ByteArrayOutputStream baos = new ByteArrayOutputStream();
    assertEquals(0, LicenseCommand.mainInit(baos, CoreCommand.INFO));

    final String result = baos.toString();
    assertTrue(result.contains("License: "));
    int longestUsageLength = "Command name".length();
    // Get longest string lengths for use below in pretty-printing.
    for (Command module : CoreCommand.INFO.commands()) {
      if (module.getCommandName().length() > longestUsageLength && LicenseCommand.showModule(module)) {
        longestUsageLength = module.getCommandName().length();
      }
    }
    final StringBuilder extraSpaces = new StringBuilder();
    for (int i = 0; i < longestUsageLength - "Command name".length(); ++i) {
      extraSpaces.append(" ");
    }
    assertTrue(result.contains("\tCommand name" + extraSpaces + " \tLicensed?         Release Level" + StringUtils.LS));
    for (Command module : CoreCommand.INFO.commands()) {
      if (module.isLicensed() || !module.isHidden()) {
        assertTrue(module.getCommandName() + " not found", contains(result, module.getCommandName(), module.isLicensed(), module.getReleaseLevel(), longestUsageLength));
        assertTrue(module.getCategory().toString() + " not found", result.contains(StringUtils.LS + StringUtils.LS + module.getCategory().getLabel() + ":" + StringUtils.LS));
      } else {
        assertFalse(module.getCommandName() + " found!", contains(result, module.getCommandName(), false, module.getReleaseLevel(), longestUsageLength));
      }
    }
  }

  private static boolean contains(String haystack, String needle, boolean licensed, ReleaseLevel releaseLevel, int longestUsageLength) {
    final StringBuilder extraSpaces = new StringBuilder();
    for (int i = 0; i < longestUsageLength - needle.length(); ++i) {
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
    assertEquals(6, LicenseCommand.getLongestLengthModule(new Command[]{ToolsCommand.FORMAT}, "abc"));
  }

}
