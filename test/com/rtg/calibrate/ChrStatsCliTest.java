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

package com.rtg.calibrate;

import java.io.File;
import java.io.IOException;

import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.AbstractCliTest;
import com.rtg.util.TestUtils;
import com.rtg.util.io.TestDirectory;

/**
 * Test class
 */
public class ChrStatsCliTest extends AbstractCliTest {

  @Override
  protected AbstractCli getCli() {
    return new ChrStatsCli();
  }

  public void testFlagValidator() throws IOException {
    try (TestDirectory td = new TestDirectory("cctest")) {
      final File template = new File(td, "in.sdf");
      final File bam = new File(td, "in.bam");
      final File ped = new File(td, "in.ped");

      checkHandleFlags("-t", template.getPath(), bam.getPath());
      checkHandleFlags("-t", template.getPath(), "--sample=foo", bam.getPath());
      checkHandleFlags("-t", template.getPath(), "-p", ped.getPath(), "--sample=foo", bam.getPath());

      final String err = checkHandleFlagsErr("-t", template.getPath(), "-p", ped.getPath(), "--sex=male", bam.getPath());
      TestUtils.containsAll(err, "Only one of --sex or --pedigree can be set");
    }
  }
}
