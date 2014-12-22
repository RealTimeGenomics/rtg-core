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
package com.rtg.variant;

import java.io.File;
import java.io.IOException;

import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.AbstractCliTest;
import com.rtg.util.StringUtils;
import com.rtg.util.io.TestDirectory;
import com.rtg.util.test.FileHelper;

/**
 *
 *
 */
public class VcfStatsCliTest extends AbstractCliTest {


  @Override
  protected AbstractCli getCli() {
    return new VcfStatsCli();
  }

  public void testInitFlags() {
    checkHelp("rtg vcfstats"
            , "Display statistics"
            , "FILE", "input VCF files"
            , "lengths", "output variant length histogram"
            );
  }

  public void testStatsRun() throws IOException {
    try (TestDirectory dir = new TestDirectory("vcfstats")) {
      final File posVcf = new File(dir, "vcfstats.vcf");
      FileHelper.resourceToFile("com/rtg/variant/resources/vcfstats.vcf", posVcf);

      checkMainInitBadFlags();
      checkMainInitWarn(posVcf.toString(), "--sample", "nosuchsample");
      String output = checkMainInitOk(posVcf.toString());
      mNano.check("vcfstats-run1.txt", StringUtils.grepMinusV(output, "Location"));

      output = checkMainInitOk(posVcf.toString(), "--sample", "sm_son1");
      mNano.check("vcfstats-run2.txt", StringUtils.grepMinusV(output, "Location"));

      output = checkMainInitOk(posVcf.toString(), "--allele-lengths");
      mNano.check("vcfstats-run3.txt", StringUtils.grepMinusV(output, "Location"));
    }

  }

}
