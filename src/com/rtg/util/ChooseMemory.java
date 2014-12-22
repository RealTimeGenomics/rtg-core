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
package com.rtg.util;

import com.rtg.util.cli.CFlags;
import com.rtg.util.cli.Validator;


/**
 * Class to choose a percentage of the available RAM for use in the RTG wrapper script.
 */
public final class ChooseMemory {

  private ChooseMemory() { }

  /**
   * Main method for getting RAM.
   * @param args the command line arguments.
   */
  public static void main(String[] args) {
    final CFlags flags = new CFlags("ChooseMemory", "Program to get the appropriate RAM to use", System.out, System.err);
    flags.registerRequired(Integer.class, "INT", "Percentage of RAM to use (1-100)");
    flags.setValidator(new Validator() {
      @Override
      public boolean isValid(CFlags flags) {
        final int percentage = (Integer) flags.getAnonymousValue(0);
        if (percentage < 1) {
          flags.setParseMessage("Percentage must be greater than 0.");
          return false;
        } else if (percentage > 100) {
          flags.setParseMessage("Percentage must be less than or equal to 100.");
          return false;
        }
        return true;
      }
    });
    if (!flags.setFlags(args)) {
      return;
    }

    final double percentage = ((Integer) flags.getAnonymousValue(0)) / 100.0;
    final int megs = (int) (Environment.getTotalMemory() * percentage / 1024.0 / 1024.0);
    System.out.println((megs < 1024 ? 1024 : megs) + "m");
  }
}
