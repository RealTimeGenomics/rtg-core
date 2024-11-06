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

package com.rtg.variant.bayes.multisample.family;

import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.IOException;

import com.rtg.launcher.AbstractNanoTest;
import com.rtg.launcher.MainResult;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.tabix.TabixIndexer;
import com.rtg.util.StringUtils;
import com.rtg.util.TestUtils;
import com.rtg.util.Utils;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.TestDirectory;
import com.rtg.util.test.BgzipFileHelper;
import com.rtg.util.test.FileHelper;
import com.rtg.variant.bayes.multisample.AbstractMultisampleCli;
import com.rtg.variant.bayes.multisample.population.PopulationNanoTest;

/**
 */
public class FamilyNanoTest extends AbstractNanoTest {

  private static final String RESOURCES_DIR = "com/rtg/variant/bayes/multisample/family/resources/";

  private static final String FAMILY_PED = ""
      + "0\tsm_dad\t0\t0\t1\t0" + StringUtils.LS
      + "0\tsm_mom\t0\t0\t2\t0" + StringUtils.LS
      + "0\tsm_son1\tsm_dad\tsm_mom\t1\t0" + StringUtils.LS
      + "0\tsm_son2\tsm_dad\tsm_mom\t1\t0" + StringUtils.LS
      + "0\tsm_daughter1\tsm_dad\tsm_mom\t2\t0" + StringUtils.LS
      + "0\tsm_daughter2\tsm_dad\tsm_mom\t2\t0" + StringUtils.LS;

  private static final String DENOVO_REF = ">simulatedSequence1 100 1" + StringUtils.LS
      + "CTTGTATCGTGTGTAATGAAAGCGCACTCATCATAGTTCTTTTGAGTCGATACTTAGCTCCCATGTGATGCGGATGACCGGTGAACC"
      + "GCAATATACATGAGTTTTAGGAACGGTACTAGTACCTCTAATACGCACTACTTTGTTCGATTACCATGGATAATCCGCCCTTGTTGG"
      + "ATTGCAATCGAACGGCAGGATCGCT" + StringUtils.LS;

  public void testNoDenovo() throws Exception {
    final String sam = FileHelper.resourceToString(RESOURCES_DIR + "family_denovo.sam");

    try (TestDirectory dir = new TestDirectory()) {
      final File priors = PopulationNanoTest.makeDenovoPriors(dir, 0.0, 0.0);
      check(DENOVO_REF, sam, FAMILY_PED, "nodenovo", new String[] {"--Xpriors", priors.getPath(), "--Xcontrary-probability", "1"});
    }
  }

  public void testDenovo() throws Exception {
    final String sam = FileHelper.resourceToString(RESOURCES_DIR + "family_denovo.sam");
    try (TestDirectory dir = new TestDirectory()) {
      final File priors = PopulationNanoTest.makeDenovoPriors(dir, 0.002, 0.002);
      check(DENOVO_REF, sam, FAMILY_PED, "denovo", new String[] {"--Xpriors", priors.getPath(), "--Xcontrary-probability", "1"});
    }
  }

  public void testDenovoHom() throws Exception {
    final String sam = FileHelper.resourceToString(RESOURCES_DIR + "family_denovo.sam");
    try (TestDirectory dir = new TestDirectory()) {
      final File priors = PopulationNanoTest.makeDenovoPriors(dir, 0.002, 0.0);
      check(DENOVO_REF, sam, FAMILY_PED, "denovo_hom", new String[] {"--Xpriors", priors.getPath(), "--Xcontrary-probability", "1"});
    }
  }

  private void check(final String template, final String samString, final String ped, final String expectedPrefix, String[] extraArgs) throws Exception {
    //Diagnostic.setLogStream(System.err);
    final File tmp = FileUtils.createTempDir("familynano", "tests");
    try {
      final File templateDir = new File(tmp, "template");
      ReaderTestUtils.getDNADir(template, templateDir);
      final File relationsFile = new File(tmp, "family.ped");
      final File out = new File(tmp, "output");
      FileUtils.stringToFile(ped, relationsFile);
      final File sam = createIndexedSamFile(samString, tmp, "total");


      final String[] baseArgs = {
          "-t", templateDir.getPath(),
          "-p", relationsFile.getPath(),
          "-o", out.getPath(),
          "-Z",
          sam.getPath(),
          "--" + AbstractMultisampleCli.NO_CALIBRATION
      };

      final String[] args;
      if (extraArgs != null) {
        args = Utils.append(baseArgs, extraArgs);
      } else {
        args = baseArgs;
      }

      final FamilyCli cli = new FamilyCli();
      final MainResult r = MainResult.run(cli, args);
      assertEquals(r.err(), 0, r.rc());

      final String result = FileUtils.fileToString(new File(out, "snps.vcf"));
      final String actualFixed = TestUtils.stripVcfHeader(result);

      //System.err.println(actualFixed);
      //final String expectedFinal = FileHelper.resourceToString("com/rtg/variant/bayes/multisample/cancer/resources/somaticnanotest" + expectedPrefix + ".vcf").replaceAll("\r", "");
      //assertEquals(expectedFinal, actualFixed.replaceAll("\r", ""));
      mNano.check("familynanotest-" + expectedPrefix + ".vcf", actualFixed, false);
      final String usageLog = cli.usageLog();
      //System.err.println(usageLog);
      TestUtils.containsAll(usageLog, "[Usage beginning module=family runId=", ", Usage end module=family runId=", " metric=9400 success=true]");
    } finally {
      assertTrue(FileHelper.deleteAll(tmp));
    }
  }

  private File createIndexedSamFile(final String cancer, final File tmp, final String tag) throws Exception {
    final File file = new File(tmp, tag + ".sam.gz");
    BgzipFileHelper.streamToBgzipFile(new ByteArrayInputStream(cancer.getBytes()), file);
    new TabixIndexer(file, new File(tmp, tag + ".sam.gz.tbi")).saveSamIndex();
    return file;
  }

  //To replicate bug 1550
  public void testBug1550() throws IOException {
    //while (true) {
    try (final TestDirectory dir = new TestDirectory()) {
      final File bam = FileHelper.resourceToFile(RESOURCES_DIR + "input.bam", new File(dir, "input.bam"));
      FileHelper.resourceToFile(RESOURCES_DIR + "input.bam.bai", new File(dir, "input.bam.bai"));
      FileHelper.resourceToFile(RESOURCES_DIR + "input.bam.calibration", new File(dir, "input.bam.calibration"));
      final File ped = FileHelper.resourceToFile(RESOURCES_DIR + "WMD.ped", new File(dir, "WMD.ped"));
      final File template = ReaderTestUtils.getDNADir(FileHelper.resourceToString(RESOURCES_DIR + "subset_template.fasta"), new File(dir, "subset_template"));
      final File outputDir = new File(dir, "output");
      final String[] args = {
          "-t", template.getPath(),
          "-o", outputDir.getPath(),
          "-p", ped.getPath(),
          bam.getPath(),
          "--region", "GL000220.1",
          "--Xignore-incompatible-sam-headers" //bams are annoying to clean up, the input contains sequence dictionary entries not in the reference :(
      };
      final MainResult r = MainResult.run(new FamilyCli(), args);
      assertEquals(r.err(), 0, r.rc());
    }
    //}
  }
}
