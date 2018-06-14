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
package com.rtg.sam;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.util.Arrays;
import java.util.Collection;

import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.AbstractCliTest;
import com.rtg.launcher.globals.GlobalFlags;
import com.rtg.util.Resources;
import com.rtg.util.SingletonPopulatorFactory;
import com.rtg.util.StringUtils;
import com.rtg.util.TestUtils;
import com.rtg.util.cli.CFlags;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.util.io.TestDirectory;
import com.rtg.util.test.FileHelper;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.util.CloseableIterator;

/**
 *
 */
public class Sam2BamTest extends AbstractCliTest {

  @Override
  protected AbstractCli getCli() {
    return new Sam2Bam();
  }

  private static final String SAMRESOURCE_1 = "com/rtg/sam/resources/mergemated.sam.gz";
  private static final String SAMRESOURCE_2 = "com/rtg/sam/resources/mergeunmated.sam.gz";
  private static final String SAMRESOURCE_3 = "com/rtg/sam/resources/mergeunmapped.sam.gz";

  /**
   * Test of convertSamToBam method, of class BamConverter.
   */
  public void testConvertSamToBam() throws Exception {
    try (final TestDirectory dir = new TestDirectory("sam2bam")) {
      final File sam1 = new File(dir, "sam1.sam.gz");
      final File sam2 = new File(dir, "sam2.sam.gz");
      final File sam3 = new File(dir, "sam3.sam.gz");
      InputStream is = Resources.getResourceAsStream(SAMRESOURCE_1);
      try {
        FileHelper.streamToFile(is, sam1);
      } finally {
        is.close();
      }
      is = Resources.getResourceAsStream(SAMRESOURCE_2);
      try {
        FileHelper.streamToFile(is, sam2);
      } finally {
        is.close();
      }
      is = Resources.getResourceAsStream(SAMRESOURCE_3);
      try {
        FileHelper.streamToFile(is, sam3);
      } finally {
        is.close();
      }
      final Collection<File> input = Arrays.asList(sam1, sam2);
      final File output = new File(dir, "output.bam");
      final File index = new File(dir, "output.bai");
      Sam2Bam.convertSamToBam(output, index, input);
      try (SamReader reader = SamUtils.makeSamReader(output)) {
        try (CloseableIterator<SAMRecord> bIt = reader.iterator()) {
          try (RecordIterator<SAMRecord> multi = new ThreadedMultifileIterator<>(input, new SingletonPopulatorFactory<>(new SamRecordPopulator()))) {
            while (bIt.hasNext() && multi.hasNext()) {
              final SAMRecord b = bIt.next();
              final SAMRecord m = multi.next();
              assertEquals(b.getSAMString().trim(), m.getSAMString().trim());
            }
            assertFalse(bIt.hasNext() || multi.hasNext());
          }
        }
      }

      GlobalFlags.resetAccessedStatus();
      final Sam2Bam bc = new Sam2Bam();
      assertTrue(output.delete());
      try (MemoryPrintStream mps = new MemoryPrintStream()) {
        final int code = bc.mainInit(new String[]{"-o", output.toString(), sam1.toString()}, mps.outputStream(), mps.printStream());
        assertEquals(mps.toString(), 0, code);
      }
      assertTrue(new File(output.toString() + ".bai").exists());
    }
  }

  public void testInnerClass() throws IOException {
    try (final TestDirectory dir = new TestDirectory("bamconverter")) {
      final String err = checkMainInitBadFlags("non-existingfile", "-o", "something");
      TestUtils.containsAll(err, "non-existingfile");

      TestUtils.containsAllUnwrapped(checkMainInitBadFlags(dir.getPath(), "-o", "something"), dir.getPath());
    }
  }

  public void testMisc() throws Exception {
    final Sam2Bam bc = new Sam2Bam();
    assertEquals("sam2bam", bc.moduleName());

    final CFlags flags = new CFlags("sdnfsdf", TestUtils.getNullPrintStream(), TestUtils.getNullPrintStream());
    Sam2Bam.initFlags(flags);
    flags.setFlags();
    assertTrue(flags.getInvalidFlagMsg().contains("Error: You must provide values for -o"));

    flags.setFlags("-o", "blah");
    assertTrue(flags.getInvalidFlagMsg().contains("Error: You must provide a value for FILE"));

    flags.setFlags("--output", "blah");
    assertTrue(flags.getInvalidFlagMsg().contains("Error: You must provide a value for FILE"));

    assertTrue(flags.getUsageString().contains("SAM/BAM format files containing mapped reads"));
    assertTrue(flags.getUsageString().contains("name for output BAM file"));
    //assertTrue(flags.getUsageString().contains("Output filename for index"));

    try (TestDirectory tmpDir = new TestDirectory("sam2bam")) {
      final File f = new File(tmpDir.getPath(), "blahdir.bam");
      assertTrue(f.mkdir());
      final File f2 = File.createTempFile("blahfile", "bam", tmpDir);
      flags.setFlags("-o", f.getPath(), f2.getPath());
      TestUtils.containsAll(flags.getParseMessage(), "already exists");
      flags.setFlags("--Xforce", "-o", f.getPath(), f2.getPath());
      TestUtils.containsAll(flags.getParseMessage(), "is a directory");
    }
  }

  public void testOutOfOrder() throws Exception {
    try (final TestDirectory tmpDir = new TestDirectory()) {
      final Sam2Bam bc = new Sam2Bam();

      final File in = new File(tmpDir, "in.sam");
      FileHelper.resourceToFile("com/rtg/sam/resources/unmated-out-of-order.sam", in);

      final MemoryPrintStream mps = new MemoryPrintStream();
      assertEquals(1, bc.mainInit(new String[]{"-o", tmpDir.getPath() + StringUtils.FS + "out.bam", in.getPath()}, TestUtils.getNullOutputStream(), mps.printStream()));
      TestUtils.containsAll(mps.toString(), "Alignments added out of order in ", "out.bam", "Sort order is coordinate.", "Offending records are at [gi0:10] and [gi0:9]");
    }
  }
  /**
   * FileHelper.resourceToFile("com/rtg/variant/resources/coverage_mated.sam.gz", mated);
   * FileHelper.resourceToFile("com/rtg/variant/resources/coverage_mated.sam.gz.calibration", new File(dir, "mated.sam.gz.calibration"));
   */
  public void testCalibrationConversion() throws IOException {
    try (final TestDirectory dir = new TestDirectory("sam2bam")) {
      final File sam = new File(dir, "samFile.sam.gz");
      final File samCal = new File(sam.getParent(), sam.getName() + ".calibration");
      final File bam = new File(dir, "bamFile.bam");
      final File bamCal = new File(bam.getParent(), bam.getName() + ".calibration");
      final File bamIndex = new File(bam.getParent(), bam.getName() + ".bai");
      FileHelper.resourceToFile("com/rtg/sam/resources/calibrated_mated.sam.gz", sam);
      FileHelper.resourceToFile("com/rtg/sam/resources/calibrated_mated.sam.gz.calibration", samCal);
      final MemoryPrintStream ps = new MemoryPrintStream();
      Diagnostic.setLogStream(ps.printStream());
      Sam2Bam.convertSamToBam(bam, bamIndex, sam, bamCal);
      assertFalse(bamCal.exists());
      assertTrue(bam.isFile());
      assertTrue(bamIndex.isFile());
      assertTrue(ps.toString(), ps.toString().contains("Number of calibration files does not match number of SAM files, will not merge calibration files."));
      ps.reset();
      assertTrue(bam.delete());
      assertTrue(bamIndex.delete());
      Sam2Bam.convertSamToBam(bam, bamIndex, sam);
      assertTrue(bam.isFile());
      assertTrue(bamIndex.isFile());
      assertTrue(bamCal.isFile());
      assertEquals(FileUtils.fileToString(samCal).replaceAll("#.*", "").replaceAll("\r", ""), FileUtils.fileToString(bamCal).replaceAll("#.*", "").replaceAll("\r", ""));
      assertFalse(ps.toString(), ps.toString().contains("Number of calibration files does not match number of SAM files, will not merge calibration files."));
    }
  }

  public void testGetBamOutputFile() {
    assertEquals("file.bam", Sam2Bam.getBamOutputFile(new File("file")).getName());
    assertEquals("file.bam", Sam2Bam.getBamOutputFile(new File("file.bam")).getName());
  }
}
