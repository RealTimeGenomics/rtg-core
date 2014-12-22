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

/**
 * Entry point for all RTG command line modules
 */
public final class Slim extends AbstractCliEntry {

  /**
   * @param args command line arguments
   */
  public static void main(String[] args) {
    new Slim().mainImpl(args);
  }

  @Override
  protected int help(String[] shiftArgs, OutputStream out, PrintStream err) {
    return CoreCommand.HELP.module().mainInit(shiftArgs, out, err);
  }

  @Override
  protected Command getSlimModule(String arg) {
    return CoreCommand.INFO.findModuleWithExpansion(arg);
  }
}
