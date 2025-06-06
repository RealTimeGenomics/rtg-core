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

import static com.rtg.util.StringUtils.LS;

import java.io.File;
import java.io.IOException;

import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.AbstractCliTest;
import com.rtg.launcher.CommonFlags;
import com.rtg.launcher.MainResult;
import com.rtg.util.StringUtils;
import com.rtg.util.TestUtils;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.SimpleArchive;
import com.rtg.util.io.TestDirectory;
import com.rtg.util.test.FileHelper;

/**
 */
public class RecalibrateCliTest extends AbstractCliTest {

  @Override
  protected AbstractCli getCli() {
    return new RecalibrateCli();
  }

  static final String EXPECTED = ", calibrate v3.0" + LS
          + "@ins:RG1\t0\t0\t0\t1" + LS
          + "@mnp:RG1\t0\t85\t8\t1\t1" + LS
          + "@nh:RG1\t0\t100" + LS
          + "@covar\treadgroup\tbasequality\tsequence\tequal\tdiff\tins\tdel" + LS
          + "RG1\t20\tsimulatedSequence1\t9889\t108\t3\t0" + LS;

  public void testSomeMethod() throws IOException {
    final File dir = FileUtils.createTempDir("recal", "test");
    try {
      final File testFile = FileHelper.resourceToFile("com/rtg/sam/resources/tinyMappings.sam.gz", new File(dir, "test.sam.gz"));
      final File templateDwa = FileHelper.resourceToFile("com/rtg/sam/resources/tinyTemplate.dwa", new File(dir, "tinyTemplate.dwa"));
      final File templateDir = new File(dir, "template");
      SimpleArchive.unpackArchive(templateDwa, templateDir);
      final MainResult r = MainResult.run(getCli(), "-t", templateDir.getPath(), testFile.getPath());
      assertEquals(r.err(), 0, r.rc());
      final File calibrationFile = new File(dir, "test.sam.gz.calibration");
      final String s = StringUtils.grepMinusV(FileUtils.fileToString(calibrationFile), "^#CL");
      final int start = s.indexOf(", calibrate ");
      assertEquals(EXPECTED, s.substring(start));
    } finally {
      assertTrue(FileHelper.deleteAll(dir));
    }
  }

  public void testPreexistCalib() throws IOException {
    final File dir = FileUtils.createTempDir("recal", "test");
    try {
      final File testFile = FileHelper.resourceToFile("com/rtg/sam/resources/tinyMappings.sam.gz", new File(dir, "test.sam.gz"));
      final File templateDwa = FileHelper.resourceToFile("com/rtg/sam/resources/tinyTemplate.dwa", new File(dir, "tinyTemplate.dwa"));
      final File templateDir = new File(dir, "template");
      SimpleArchive.unpackArchive(templateDwa, templateDir);
      final File calib = new File(dir, testFile.getName() + CommonFlags.RECALIBRATE_EXTENSION);
      assertTrue(calib.createNewFile());
      MainResult r = MainResult.run(getCli(), "-t", templateDir.getPath(), testFile.getPath());
      assertEquals(1, r.rc());
      TestUtils.containsAll(r.err(), "Error: Calibration file already exists:", calib.getPath());

      r = MainResult.run(getCli(), "-t", templateDir.getPath(), testFile.getPath(), "--force");
      assertEquals(r.err(), 0, r.rc());
    } finally {
      assertTrue(FileHelper.deleteAll(dir));
    }
  }

  private static final String[] EXPECTED2 = {
    ", calibrate v3.0" + LS,
    "@ins:RG1\t0\t0\t0\t1" + LS,
    "@mnp:RG1\t0\t85\t8\t1\t1" + LS,
    "@nh:RG1\t0\t100" + LS,
    "@covar\treadgroup\tbasequality\tmachinecycle:100\tequal\tdiff\tins\tdel" + LS,
    "RG1\t20\t0\t99\t1\t0\t0" + LS,
    "RG1\t20\t1\t100\t0\t0\t0" + LS,
    "RG1\t20\t2\t98\t2\t0\t0" + LS,
    "RG1\t20\t3\t100\t0\t0\t0" + LS,
    "RG1\t20\t4\t98\t1\t1\t0" + LS,
    // ...
    "RG1\t20\t98\t99\t1\t0\t0" + LS,
    "RG1\t20\t99\t100\t0\t0\t0" + LS,
  };

  public void testCovariates() throws IOException {
    final File dir = FileUtils.createTempDir("recal", "test");
    try {
      final File testFile = FileHelper.resourceToFile("com/rtg/sam/resources/tinyMappings.sam.gz", new File(dir, "test.sam.gz"));
      final File templateDwa = FileHelper.resourceToFile("com/rtg/sam/resources/tinyTemplate.dwa", new File(dir, "tinyTemplate.dwa"));
      final File templateDir = new File(dir, "template");
      SimpleArchive.unpackArchive(templateDwa, templateDir);
      final MainResult r = MainResult.run(getCli(), "--force", "--template", templateDir.getPath(), testFile.getPath(),
          "--Xcovariate=machinecycle",
          "--Xcovariate=readgroup",
          "-c", "basequality");
      assertEquals(r.err(), 0, r.rc());
      final File calibrationFile = new File(dir, "test.sam.gz.calibration");
      final String s = StringUtils.grepMinusV(FileUtils.fileToString(calibrationFile), "^#CL");
      final int start = s.indexOf(", calibrate ");
      TestUtils.containsAll(s.substring(start), EXPECTED2);
    } finally {
      assertTrue(FileHelper.deleteAll(dir));
    }
  }
  public void testInitParams() {
    checkHelp("rtg calibrate [OPTION]... -t SDF FILE+",
        " [OPTION]... -t SDF -I FILE",
        "Creates quality calibration files for all supplied SAM/BAM files.",
        "file containing a list of SAM/BAM format files (1 per line) containing mapped reads",
        "SDF containing the reference",
        "SAM/BAM format files containing mapped reads. May be specified 0 or more times",
        "print help on command-line flag usage"
    );
  }
  private static final String EXPECTED_WITH_BED = ", calibrate v3.0" + LS
                                         + "@nh:RG1\t0\t4" + LS
                                         + "@sequence\t7\tsimulatedSequence1" + LS
                                         + "@covar\treadgroup\tbasequality\tsequence\tequal\tdiff\tins\tdel" + LS
                                         + "RG1\t20\tsimulatedSequence1\t28\t0\t0\t0" + LS;
  public void testWithBed() throws IOException {
    try (TestDirectory dir = new TestDirectory()) {
      final File testFile = FileHelper.resourceToFile("com/rtg/sam/resources/tinyMappings.sam.gz", new File(dir, "test.sam.gz"));
      final File templateDwa = FileHelper.resourceToFile("com/rtg/sam/resources/tinyTemplate.dwa", new File(dir, "tinyTemplate.dwa"));
      final File bedFile = FileHelper.resourceToFile("com/rtg/sam/resources/calibrate.bed", new File(dir, "calibrate.bed"));
      final File templateDir = new File(dir, "template");
      SimpleArchive.unpackArchive(templateDwa, templateDir);
      final MainResult r = MainResult.run(getCli(), "-t", templateDir.getPath(), testFile.getPath(), "--bed-regions", bedFile.getPath());
      assertEquals(r.err(), 0, r.rc());
      final File calibrationFile = new File(dir, "test.sam.gz.calibration");
      final String s = StringUtils.grepMinusV(FileUtils.fileToString(calibrationFile), "^#CL");
      final int start = s.indexOf(", calibrate ");
      assertEquals(EXPECTED_WITH_BED, s.substring(start));
    }
  }
}
