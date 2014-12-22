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

import com.rtg.util.Constants;
import com.rtg.util.StringUtils;
import com.rtg.util.cli.CFlags;
import com.rtg.util.cli.WrappingStringBuilder;

/**
 * Help module for Slim, provides summaries of other modules
 *
 */
public final class HelpCommand {

  private HelpCommand() { } //prevent instantiation

  private static final String APPLICATION_NAME = Constants.APPLICATION_NAME;

  private static final ResourceBundle USAGE_RESOURCE = ResourceBundle.getBundle("com.rtg.Usage", Locale.getDefault());
  static final String USAGE_STR = "Usage: " + APPLICATION_NAME + " COMMAND [OPTION]..." + StringUtils.LS + StringUtils.LS
                                + "Type '" + APPLICATION_NAME + " help COMMAND' for help on a specific command." + StringUtils.LS
                                + "The following commands are available:" + StringUtils.LS;

  /**
   * Main function, entry-point for help.
   * @param args currently unused, will take module name in future
   * @param outStream print stream for output
   * @param errStream print stream for error
   * @param info source of module command names
   * @return shell return code 0 for success, anything else for failure
   */
  public static int mainInit(final String[] args, final OutputStream outStream, final PrintStream errStream, CommandLookup info) {
    //final PrintStream outPrintStream = new PrintStream(outStream);
    return printHelp(args, outStream, errStream, info);
  }

  private static int printHelp(final String[] args, final OutputStream outStream, final PrintStream errStream, CommandLookup info) {
    if (args != null && args.length > 0) {
      final String args0 = args[0].toUpperCase(Locale.getDefault());
      if (args0.equals("--XHIDDEN") || args0.equals("--XHELP")) {
        final PrintStream psOutStream =  new PrintStream(outStream);
        try {
          printUsage(psOutStream, true, info);
        } finally {
          psOutStream.flush();
        }
        return 0;
      } else {
        final Command mod = info.findModule(args0);
        if (mod == null) {
          final PrintStream psOutStream =  new PrintStream(outStream);
          try {
            printUsage(psOutStream, false, info);
          } finally {
            psOutStream.flush();
          }
        } else {
          if (!mod.isLicensed()) {
            outputUnlicensedModule(mod);
          } else {
            mod.mainInit(new String[] {"-h"}, outStream, errStream);
          }
        }
        return 0; //successfully determined help for this module
      }
    } else {  //args
      final PrintStream psOutStream = new PrintStream(outStream);
      try {
        printUsage(args == null ? errStream : psOutStream, false, info);
      } finally {
        psOutStream.flush();
      }
      return args == null ? 1 : 0; //if no args given, return error
    }
  }

  protected static void outputUnlicensedModule(Command mod) {
    final String message = "The " + mod.getCommandName()
      + " command has not been enabled by your current license." + StringUtils.LS + "Please contact "
      + Constants.SUPPORT_EMAIL_ADDR + " to have this command licensed.";
    System.err.println(message);
    throw new RuntimeException();
  }

  /**
   * Usage information for the SLIM runtime. Includes information about modules.
   * @param printHidden print hidden modules
   * @param info source of module command names
   * @return string representation of usage information
   */
  public static String getUsage(boolean printHidden, CommandLookup info) {
    return getUsage(printHidden, CFlags.DEFAULT_WIDTH, info);
  }

  protected static String getUsage(boolean printHidden, int width, CommandLookup info) {
    String moduleTypeName;
    moduleTypeName = null;

    // Get longest string lengths for use below in pretty-printing.
    final int longestUsageLength = Math.min(getLongestLengthModule(info.commands()), 13);

    final StringBuilder sb = new StringBuilder();

    sb.append(USAGE_STR).append(StringUtils.LS);
    for (Command module : info.commands()) {

      // Show only licensed, non-hidden modules.
      if (!module.isLicensed()) {
        continue;
      }
      if (!printHidden && module.isHidden())  {
        continue;
      }

      if (moduleTypeName == null || !moduleTypeName.equals(module.getCategory().toString())) {
        if (moduleTypeName != null) {
          sb.append(StringUtils.LS);
        }
        moduleTypeName = module.getCategory().toString();
        try {
          sb.append(USAGE_RESOURCE.getString(moduleTypeName + "_CAT_DESC")).append(":").append(StringUtils.LS);
        } catch (MissingResourceException mre) {
          sb.append(moduleTypeName).append(":").append(StringUtils.LS);
        }
      }
      sb.append("\t").append(module.getCommandName().toLowerCase(Locale.getDefault()));
      if (module.getCommandName().length() > longestUsageLength) {
        sb.append(StringUtils.LS).append("\t");
        for (int i = 0; i < longestUsageLength; i++) {
          sb.append(" ");
        }
      } else {
        for (int i = 0; i < longestUsageLength - module.getCommandName().length(); i++) {
          sb.append(" ");
        }
      }
      sb.append(" \t");
      try {
        sb.append(USAGE_RESOURCE.getString(module.getCommandName() + "_CMD_DESC")).append(StringUtils.LS);
      } catch (MissingResourceException mre) {
        sb.append(StringUtils.LS);  //just a blank...
      }
    }
    final WrappingStringBuilder wb = new WrappingStringBuilder();
    wb.setWrapWidth(width);
    final StringBuilder spaces = new StringBuilder();
    for (int i = 0; i < longestUsageLength; i++) {
      spaces.append(" ");
    }
    wb.setWrapIndent("\t" + spaces.toString() + " \t");
    wb.wrapTextWithNewLines(sb.toString());
    return wb.toString();
  }

  /**
   * Print the usage information for the SLIM runtime to the given print stream
   * @param printStream the print stream to print usage information to
   * @param printHidden print hidden modules
   * @param info source of module command names
   */
  public static void printUsage(final PrintStream printStream, final boolean printHidden, CommandLookup info) {
    printStream.print(getUsage(printHidden, info));
  }

  static int getLongestLengthModule(Command[] values) {
    int longestUsageLength = 0;
    for (Command module : values) {
      if (module.getCommandName().length() > longestUsageLength && !module.isHidden()) {
        longestUsageLength = module.getCommandName().length();
      }
    }
    return longestUsageLength;
  }
}
