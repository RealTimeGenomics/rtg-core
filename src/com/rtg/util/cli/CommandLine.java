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
package com.rtg.util.cli;

import java.util.Arrays;
import java.util.UUID;

import com.rtg.util.StringUtils;

/**
 * Utility methods for storing and retrieving the command line arguments.
 *
 */
public final class CommandLine {

  private static String[] sCommandArgs;
  private static String sCommandLine;

  private static UUID sRunId = UUID.randomUUID();

  private CommandLine() {
  }

  /**
   * Sets the current command line
   *
   * @param args the command line
   */
  public static void setCommandArgs(String... args) {
    sCommandArgs = Arrays.copyOf(args, args.length);
    sCommandLine = StringUtils.implode(args, " ", true);
  }

  /**
   * Revert to default state
   */
  public static void clearCommandArgs() {
    sCommandArgs = null;
    sCommandLine = null;
  }

  /**
   * Get the command line arguments used to invoke this instance of Slim
   * @return the command line as individual arguments
   */
  public static String[] getCommandArgs() {
    return Arrays.copyOf(sCommandArgs, sCommandArgs.length);
  }

  /**
   * Get a single printable string containing the command line used to invoke this instance of Slim
   * @return the command line as a single string.
   */
  public static String getCommandLine() {
    return sCommandLine;
  }

  /**
   * @return a UUID for this instance of Slim
   */
  public static UUID getRunId() {
    return sRunId;
  }
}
