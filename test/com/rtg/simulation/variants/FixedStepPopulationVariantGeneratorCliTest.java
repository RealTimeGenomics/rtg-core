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
package com.rtg.simulation.variants;

import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.AbstractCliTest;

/**
 */
public class FixedStepPopulationVariantGeneratorCliTest extends AbstractCliTest {

  @Override
  protected AbstractCli getCli() {
    return new FixedStepPopulationVariantGeneratorCli();
  }

  /** Test for an error that will be picked up during flags parsing. */
  private void checkFlagsError(final String[] args, final String exp) {
    final String err = checkHandleFlagsErr(args);
    assertTrue("<" + exp + "> was not in <" + err + ">", err.contains(exp));
  }

  private static final String EXP_F1 = "Error: You must provide values for -d INT -o FILE -i SDF --spec STRING";

  public void testErrorFlags() {
    checkFlagsError(new String[] {}, EXP_F1);
  }


  public void testInitParams() {
    checkHelp("-d INT -o FILE -i SDF --spec STRING",
        "distance between mutations",
        "SDF containing input genome",
        "allele frequency"
    );
  }
}
