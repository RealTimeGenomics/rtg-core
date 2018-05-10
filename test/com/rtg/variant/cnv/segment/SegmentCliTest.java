/*
 * Copyright (c) 2016. Real Time Genomics Limited.
 *
 * Use of this source code is bound by the Real Time Genomics Limited Software Licence Agreement
 * for Academic Non-commercial Research Purposes only.
 *
 * If you did not receive a license accompanying this file, a copy must first be obtained by email
 * from support@realtimegenomics.com.  On downloading, using and/or continuing to use this source
 * code you accept the terms of that license agreement and any amendments to those terms that may
 * be made from time to time by Real Time Genomics Limited.
 */

package com.rtg.variant.cnv.segment;

import java.io.File;
import java.io.IOException;

import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.AbstractCliTest;
import com.rtg.launcher.MainResult;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.reader.SdfId;
import com.rtg.util.TestUtils;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.TestDirectory;
import com.rtg.util.test.FileHelper;

import htsjdk.samtools.util.BlockCompressedInputStream;

/**
 * Tests the corresponding class.
 */
public class SegmentCliTest extends AbstractCliTest {

  @Override
  protected AbstractCli getCli() {
    return new SegmentCli();
  }

  public void testHelp() {
    checkHelp("Segments depth of coverage data to identify regions of consistent copy number.",
      "BED file supplying per-region coverage data for the sample",
      "BED file supplying per-region coverage data for control sample",
      "directory for output",
      "reference genome",
      "weighting factor for inter-segment distances during energy scoring",
      "weighting factor for intra-segment distances during energy scoring",
      "segmentation sensitivity",
      "Sensitivity Tuning",
      "Reporting",
      "Utility",
      "File Input/Output",
      "--no-gzip",
      "--help"
      );
  }

  public void testErrors() throws IOException {
    final File emptyFile = File.createTempFile("test", ".vcf");
    try {
      String res = checkHandleFlagsErr();
      assertTrue(res.contains("Error: You must provide values for --case FILE -o DIR -t SDF"));
      res = checkHandleFlagsErr("-o", "test-foo-out", "-t", "test-sdf", "--case", emptyFile.getPath());
      TestUtils.containsAll(res, "Error: One of --Xcolumn or --control or --panel must be set");
      res = checkHandleFlagsErr("-o", "test-foo-out", "-t", "test-sdf", "--case", emptyFile.getPath(), "--control", emptyFile.getPath(), "--Xmin-segments", "0");
      assertTrue(res.contains("Error: The value for --Xmin-segments must be at least 1"));
    } finally {
      assertTrue(emptyFile.delete());
    }
  }

  public void testDesc() {
    assertEquals("segment depth of coverage data to identify copy number alterations", getCli().description());
  }

  public void testSimpleOperation() throws IOException {
    try (final TestDirectory dir = new TestDirectory()) {
      final File reference = new File(dir, "reference");
      ReaderTestUtils.getReaderDNA(mNano.loadReference("ref.fa"), reference, new SdfId(0));

      final File control = FileHelper.resourceToFile("com/rtg/variant/cnv/segment/resources/control.bed", new File(dir, "control.bed"));
      final File sample = FileHelper.resourceToFile("com/rtg/variant/cnv/segment/resources/case.bed", new File(dir, "case.bed"));
      final File output = new File(dir, "output");

      // gc correction does not work well here because the regions are small, so turn it off
      final MainResult result = MainResult.run(getCli(), "-t", reference.getPath(), "-o", output.getPath(), "--control", control.getPath(), "--case", sample.getPath(), "--sample", "foo", "--beta", "0.1", "--Xgcbins", "0");
      assertEquals(0, result.rc());
      final File vcf = new File(output, "segments.vcf.gz");
      assertEquals(BlockCompressedInputStream.FileTermination.HAS_TERMINATOR_BLOCK, BlockCompressedInputStream.checkTermination(vcf));
      mNano.check("expected.unsegmented.bed", FileHelper.gzFileToString(new File(output, "unsegmented.bed.gz")));
      mNano.check("expected.segments.vcf", TestUtils.sanitizeVcfHeader(FileHelper.gzFileToString(vcf)));
      mNano.check("expected.summary.txt", FileUtils.fileToString(new File(output, "summary.txt")));
      mNano.check("expected.out.txt", result.out());
    }
  }

  public void testPanelOperation() throws IOException {
    try (final TestDirectory dir = new TestDirectory()) {
      final File reference = new File(dir, "reference");
      ReaderTestUtils.getReaderDNA(mNano.loadReference("ref.fa"), reference, new SdfId(0));

      final File panel = FileHelper.resourceToFile("com/rtg/variant/cnv/segment/resources/panel.bed", new File(dir, "panel.bed"));
      final File sample = FileHelper.resourceToFile("com/rtg/variant/cnv/segment/resources/case.bed", new File(dir, "case.bed"));
      final File output = new File(dir, "output");

      // gc correction does not work well here because the regions are small, so turn it off
      final MainResult result = MainResult.run(getCli(), "-t", reference.getPath(), "-o", output.getPath(), "--panel", panel.getPath(), "--case", sample.getPath(), "--sample", "foo", "-Z", "--beta", "0.1", "--Xgcbins", "0", "--min-control-coverage", "0.001");
      assertEquals(0, result.rc());
      mNano.check("expected.unsegmented-panel.bed", FileUtils.fileToString(new File(output, "unsegmented.bed")));
      mNano.check("expected.segments-panel.vcf", TestUtils.sanitizeVcfHeader(FileUtils.fileToString(new File(output, "segments.vcf"))));
      mNano.check("expected.summary-panel.txt", FileUtils.fileToString(new File(output, "summary.txt")));
      mNano.check("expected.out-panel.txt", result.out());
    }
  }

  public void testGcBiasOperation() throws IOException {
    // This is not the best test because the bins available in this small test are rather
    // unpopulated, but it does at least exercise some of the GC code.
    try (final TestDirectory dir = new TestDirectory()) {
      final File reference = new File(dir, "reference");
      ReaderTestUtils.getReaderDNA(mNano.loadReference("ref.fa"), reference, new SdfId(0));

      final File control = FileHelper.resourceToFile("com/rtg/variant/cnv/segment/resources/control.bed", new File(dir, "control.bed"));
      final File sample = FileHelper.resourceToFile("com/rtg/variant/cnv/segment/resources/case.bed", new File(dir, "case.bed"));
      final File output = new File(dir, "output");

      // gc correction does not work well here because the regions are small, so turn it off
      final MainResult result = MainResult.run(getCli(), "-t", reference.getPath(), "-o", output.getPath(), "--control", control.getPath(), "--case", sample.getPath(), "--sample", "foo", "-Z", "--beta", "0.1", "--Xgcbins", "2");
      assertEquals(0, result.rc());
      mNano.check("expected.unsegmented-gc.bed", FileUtils.fileToString(new File(output, "unsegmented.bed")));
      mNano.check("expected.segments-gc.vcf", TestUtils.sanitizeVcfHeader(FileUtils.fileToString(new File(output, "segments.vcf"))));
      mNano.check("expected.summary-gc.txt", FileUtils.fileToString(new File(output, "summary.txt")));
      mNano.check("expected.out-gc.txt", result.out());
    }
  }

}
