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

import static com.rtg.util.StringUtils.LS;

import java.io.File;
import java.io.IOException;

import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.AbstractCliTest;
import com.rtg.util.TestUtils;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.MemoryPrintStream;
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
          + "#CL\tnull" + LS
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
      final MemoryPrintStream dump = new MemoryPrintStream();
      final int code = getCli().mainInit(new String[] {"-t", templateDir.getPath(), testFile.getPath()}, dump.outputStream(), dump.printStream());
      assertEquals(dump.toString(), 0, code);
      final File calibrationFile = new File(dir, "test.sam.gz.calibration");
      final String s = FileUtils.fileToString(calibrationFile);
      System.err.println(s);
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
      final File calib = new File(dir, testFile.getName() + Recalibrate.EXTENSION);
      assertTrue(calib.createNewFile());
      final MemoryPrintStream dump = new MemoryPrintStream();
      int code = getCli().mainInit(new String[] {"-t", templateDir.getPath(), testFile.getPath()}, dump.outputStream(), dump.printStream());
      assertEquals(1, code);
      TestUtils.containsAll(dump.toString(), "Error: Calibration file already exists:", calib.getPath());

      code = getCli().mainInit(new String[] {"-t", templateDir.getPath(), testFile.getPath(), "--force"}, dump.outputStream(), dump.printStream());
      assertEquals(0, code);
    } finally {
      assertTrue(FileHelper.deleteAll(dir));
    }
  }

  private static final String[] EXPECTED2 = {
    ", calibrate v3.0" + LS,
    "#CL\tnull" + LS,
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
      final MemoryPrintStream dump = new MemoryPrintStream();
      final int code = getCli().mainInit(new String[] {"--force", "--template", templateDir.getPath(), testFile.getPath(),
          "--Xcovariate=machinecycle",
          "--Xcovariate=readgroup",
          "-c", "basequality"},
          dump.outputStream(), dump.printStream());
      assertEquals(dump.toString(), 0, code);
      final File calibrationFile = new File(dir, "test.sam.gz.calibration");
      final String s = FileUtils.fileToString(calibrationFile);
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
        "SDF containing reference",
        "SAM/BAM format files containing mapped reads. May be specified 0 or more times",
        "force overwriting of calibration files",
        "print help on command-line flag usage"
    );
  }
  private static final String EXPECTED_WITH_BED = ", calibrate v3.0" + LS
                                         + "#CL\tnull" + LS
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
      try (MemoryPrintStream dump = new MemoryPrintStream()) {
        final int code = getCli().mainInit(new String[]{"-t", templateDir.getPath(), testFile.getPath(), "--bed-regions", bedFile.getPath()}, dump.outputStream(), dump.printStream());
        assertEquals(dump.toString(), 0, code);
        final File calibrationFile = new File(dir, "test.sam.gz.calibration");
        final String s = FileUtils.fileToString(calibrationFile);
        final int start = s.indexOf(", calibrate ");
        assertEquals(EXPECTED_WITH_BED, s.substring(start));
      }
    }
  }
}
