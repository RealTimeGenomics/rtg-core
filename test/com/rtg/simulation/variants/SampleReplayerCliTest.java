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

import java.io.File;
import java.io.IOException;

import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.AbstractCliTest;

/**
 * Test the corresponding class
 */
public class SampleReplayerCliTest extends AbstractCliTest {

  protected File mDir;

  @Override
  public void setUp() throws IOException {
    super.setUp();
  }

  @Override
  public void tearDown() throws IOException {
    super.tearDown();
  }

  @Override
  protected AbstractCli getCli() {
    return new SampleReplayerCli();
  }

  public void testInitParams() {
    checkHelp("reference", "reference genome",
        "-i", "input=", "input VCF",
        "-o", "output=", "output SDF",
        "print help on command-line flag usage",
        "sample=", "name of the sample");
  }

}
