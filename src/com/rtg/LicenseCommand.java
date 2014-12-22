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

import java.io.OutputStream;
import java.io.PrintStream;
import java.util.Locale;
import java.util.MissingResourceException;
import java.util.ResourceBundle;

import com.rtg.util.License;
import com.rtg.util.StringUtils;

/**
 * License command for Slim, provides summaries of other commands
 *
 */
public final class LicenseCommand {

  private static final ResourceBundle USAGE_RESOURCE = ResourceBundle.getBundle("com.rtg.Usage", Locale.getDefault());

  private LicenseCommand() { } //prevent instantiation

  /**
   * Main function, entry-point for license.
   *
   * @param outStream print stream for output
   * @param info source of module command names
   * @return shell return code 0 for success, anything else for failure
   */
  public static int mainInit(final OutputStream outStream, CommandLookup info) {
    final PrintStream psoutStream = new PrintStream(outStream);
    try {
      printLicense(psoutStream, info);
    } finally {
      psoutStream.flush();
    }
    return 0;
  }

  static String getLicenseSummary() {
    final StringBuilder sb = new StringBuilder();
    sb.append("License: ").append(License.getMessage()).append(StringUtils.LS);
    final String person = License.getPerson();
    if (person != null) {
      sb.append("Licensed to: ").append(person).append(StringUtils.LS);
    }
    final String keyPath = License.getKeyPath();
    if (keyPath != null) {
      sb.append("License location: ").append(keyPath).append(StringUtils.LS);
    }
    return sb.toString();
  }

  private static String sModuleTypeName;

  static void resetModuleTypeName() {
    sModuleTypeName = null;
  }

  private static void printLicense(final PrintStream out, CommandLookup info) {
    resetModuleTypeName();

    out.println(getLicenseSummary());

    final String commandName = "Command name";
    // Get longest string lengths for use below in pretty-printing.
    final int longestUsageLength = getLongestLengthModule(info.commands(), commandName);
    out.print("\t");
    out.print(commandName);
    for (int i = -1; i < longestUsageLength - commandName.length(); i++) {
      out.print(" ");
    }
    out.print("\t" + padTo("Licensed?", 17));
    out.print(" Release Level");
    out.println();
    out.println();

    for (Command module : info.commands()) {
      outputModule(module, longestUsageLength, out);
    }
  }

  static void outputModule(Command module, int longestUsageLength, PrintStream out) {
    if (showModule(module)) {
      if (sModuleTypeName == null || !sModuleTypeName.equals(module.getCategory().toString())) {
        if (sModuleTypeName != null) {
          out.println();
        }
        sModuleTypeName = module.getCategory().toString();
        try {
          out.println(USAGE_RESOURCE.getString(sModuleTypeName + "_CAT_DESC") + ":");
        } catch (MissingResourceException mre) {
          out.println(sModuleTypeName + ":");
        }
      }
      out.print("\t" + module.getCommandName().toLowerCase(Locale.getDefault()));
      for (int i = 0; i < longestUsageLength - module.getCommandName().length(); i++) {
        out.print(" ");
      }

      out.print(" \t");
      out.print(padTo(module.licenseStatus(), 17));

      out.print(" ");
      if (module.getReleaseLevel() == ReleaseLevel.GA) {
        out.print(module.getReleaseLevel());
      } else {
        out.print(module.getReleaseLevel().toString().toLowerCase(Locale.getDefault()));
      }

      out.println();
    }
  }

  static String padTo(String str, int length) {
    final StringBuilder sb = new StringBuilder(str);
    for (int i = str.length(); i < length; i++) {
      sb.append(' ');
    }
    return sb.toString();
  }

  // Show all GA and beta modules modules as they have had some level of polish applied.
  // Developers also get to see warty alpha modules
  static boolean showModule(final Command module) {
    return (module.getReleaseLevel() != ReleaseLevel.ALPHA) || License.isDeveloper();
  }

  static int getLongestLengthModule(Command[] values, String baseString) {
    int longestUsageLength = baseString.length();
    for (Command module : values) {
      if (module.getCommandName().length() > longestUsageLength && showModule(module)) {
        longestUsageLength = module.getCommandName().length();
      }
    }
    return longestUsageLength;
  }
}
