/*
 * Copyright (c) 2018. Real Time Genomics Limited.
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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

  @Override
  public void testDescription() {
    assertEquals(ChrStatsCli.DESCRIPTION, getCli().description());
  }

  public void testInitParams() {
    checkHelp("chrstats [OPTION]... -t SDF FILE+"
      , "-t SDF -I FILE"
      , "Check expected chromosome coverage levels from mapping calibration files."
      , "File Input/Output"
      , "-t,", "--template=SDF", "reference genome"
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
