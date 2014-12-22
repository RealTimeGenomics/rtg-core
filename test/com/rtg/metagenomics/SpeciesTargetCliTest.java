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
package com.rtg.metagenomics;

import java.io.File;
import java.io.IOException;

import com.rtg.launcher.AbstractParamsCliTest;
import com.rtg.launcher.GlobalFlags;
import com.rtg.launcher.ParamsCli;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.reader.SdfId;
import com.rtg.util.TestUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.FileUtils;
import com.rtg.util.test.FileHelper;

/**
 */
public class SpeciesTargetCliTest extends AbstractParamsCliTest<SpeciesParams> {

  private static final String APP_NAME = "rtg species-target";

  @Override
  public void setUp() throws IOException {
    super.setUp();
    GlobalFlags.resetAccessedStatus();
    //comment this out if you want diagnostics from the test
    Diagnostic.setLogStream();
  }

  @Override
  public final void testApplicationName() {
    assertEquals(APP_NAME, new SpeciesTargetCli().applicationName() + " " + new SpeciesTargetCli().moduleName());
  }

  public final void testInitFlags() {
    checkHelp(APP_NAME
        , "-t,", "--genomes=SDF", "SDF containing the genomes"
        , "-I,", "--input-list-file=FILE", "file containing a list of SAM/BAM format files (1 per line) containing mapped reads"
        , "-o,", "--output=DIR", "directory for output"
        , "--exclude-mated", "exclude all mated SAM records"
        , "--exclude-unmated", "exclude all unmated SAM records"
        , "-h,", "--help", "print help on command-line flag usage"
        );
  }

  @Override
  protected ParamsCli<SpeciesParams> getParamsCli() {
    return new SpeciesTargetCli();
  }

  public void testEndToEnd() throws IOException {
    final File dir = FileUtils.createTempDir("test", "speciesCli");
    try {
      final File t = new File(dir, "t");
      ReaderTestUtils.getReaderDNA(">a1\nacgt\n>a2\nacgt\n>b1\nacgt\n>b2\nacgt\n", t, new SdfId(1L)).close();
      final File o = new File(dir, "o");
      final File r = new File(dir, "r");
      FileUtils.stringToFile("a1\ta n a\na2\ta n a\nb1\tb\nb2\tb", r);
      final File s = new File(dir, "s.sam");
      FileUtils.stringToFile(SpeciesCliTest.EXTRA_SAM, s);
      checkMainInitOk("-o", o.getPath(), "-t", t.getPath(), s.getPath(), "--target", "0");
      final String usage = mCli.usageLog();
      TestUtils.containsAll(usage, "Usage beginning module=species-target runId=", ", Usage end module=species-target runId=", " metric=6 success=true");
    } finally {
      assertTrue(FileHelper.deleteAll(dir));
    }
  }

}

