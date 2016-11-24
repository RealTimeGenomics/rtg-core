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

import static org.hamcrest.core.StringStartsWith.startsWith;
import static org.junit.Assert.assertThat;

import java.io.File;
import java.io.IOException;

import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.AbstractCliTest;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.util.TestUtils;
import com.rtg.util.io.FileUtils;
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
      TestUtils.containsAll(err, "Cannot set both --sex and --pedigree");
      if (!bam.createNewFile() && ped.createNewFile()) {
        throw new IOException();
      }
      ReaderTestUtils.getDNADir(">simulatedSequence1\nacgtacgtacgtacgtacgt", template);

      assertThat(checkMainInitBadFlags("-t", template.getPath(), "-p", ped.getPath(), "--sample=foo", bam.getPath()),
        startsWith("Error: Number of calibration files (0) does not match number of SAM files (1)"));

      FileUtils.stringToFile("#Version 1" + RecalibrateCliTest.EXPECTED, new File(td, "in.bam.calibration"));
      final String s = checkMainInitBadFlags("-t", template.getPath(), "-p", ped.getPath(), "--sample=foo", bam.getPath());
      assertThat(s, startsWith("Error: The supplied reference does not contain sufficient genome chromosome information."));
    }
  }

  public void testDescription() {
    assertEquals(ChrStatsCli.DESCRIPTION, getCli().description());
  }

  public void testInitParams() {
    checkHelp("chrstats [OPTION]... -t SDF FILE+"
      , "-t SDF -I FILE"
      , "Check expected chromosome coverage levels from mapping calibration files."
      , "File Input/Output"
      , "-t,", "--template=SDF", "SDF containing reference genome"
      , "-I,", "--input-list-file=FILE", "file containing a list of SAM/BAM format files (1 per line) containing mapped reads"
      , "-p,", "--pedigree=FILE", "genome relationships PED file"
      , "--sex=SEX", "sex setting that the individual was mapped as"
      , "-s,", "--sample=STRING", "the name of the sample to check"
    );
    checkExtendedHelp("chrstats [OPTION]... -t SDF FILE+"
        , "-t SDF -I FILE"
      , "--Xhelp", "print help on extended command-line flag usage"
    );
  }
}
