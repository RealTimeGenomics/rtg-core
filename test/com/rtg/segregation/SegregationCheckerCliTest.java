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

package com.rtg.segregation;

import java.io.File;
import java.io.IOException;

import com.rtg.reader.ReaderTestUtils;
import com.rtg.reference.ReferenceGenome;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.Talkback;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.util.io.TestDirectory;
import com.rtg.util.test.FileHelper;
import com.rtg.util.test.NanoRegression;

import junit.framework.TestCase;

/**
 */
public class SegregationCheckerCliTest extends TestCase {

  private NanoRegression mNano = null;

  @Override
  public void setUp() {
    Diagnostic.setLogStream();
    mNano = new NanoRegression(this.getClass(), false);
  }

  @Override
  public void tearDown() throws Exception {
    // clear the module name so later tests don't report SlimException to the
    // Talkback system
    Talkback.setModuleName(null);
    try {
      mNano.finish();
    } finally {
      mNano = null;
    }
  }

  public void testNano() throws IOException {
    try (final TestDirectory dir = new TestDirectory()) {
      final File template = ReaderTestUtils.getDNADir(">1\nACGT\n>2\nACGT\n>3\nACGT", new File(dir, "template"));
      FileUtils.stringToFile("version 1\neither def diploid linear", new File(template, ReferenceGenome.REFERENCE_FILE));
      final MemoryPrintStream out = new MemoryPrintStream();
      final MemoryPrintStream err = new MemoryPrintStream();
      final File vcf = FileHelper.resourceToFile("com/rtg/segregation/resources/crossover.vcf", new File(dir, "vcf.vcf"));
      final File bed = FileHelper.resourceToFile("com/rtg/segregation/resources/regions.bed", new File(dir, "regions.bed"));
      final File output = new File(dir, "out.vcf");
      final SegregationCheckerCli cli = new SegregationCheckerCli();
      final String[] args = {"--template", template.getPath(), "--bed", bed.getPath(), "--vcf", vcf.getPath(), "--output", output.getPath(), "--father", "Father", "--mother", "Mother", "-Z"};
      final int ret = cli.mainInit(args, out.outputStream(), err.printStream());
      assertEquals(err.toString(), 0, ret);
      final String result = FileUtils.fileToString(output);
      mNano.check("crossovers_annotated.vcf", result, true);
    }
  }

  public void testNanoRepair() throws IOException {
    try (final TestDirectory dir = new TestDirectory()) {
      final File template = ReaderTestUtils.getDNADir(">1\nACGT\n>2\nACGT\n>3\nACGT", new File(dir, "template"));
      FileUtils.stringToFile("version 1\neither def diploid linear", new File(template, ReferenceGenome.REFERENCE_FILE));
      final MemoryPrintStream out = new MemoryPrintStream();
      final MemoryPrintStream err = new MemoryPrintStream();
      final File vcf = FileHelper.resourceToFile("com/rtg/segregation/resources/crossover.vcf", new File(dir, "vcf.vcf"));
      final File bed = FileHelper.resourceToFile("com/rtg/segregation/resources/regions.bed", new File(dir, "regions.bed"));
      final File output = new File(dir, "out.vcf");
      final SegregationCheckerCli cli = new SegregationCheckerCli();
      final String[] args = {"--repair", "--template", template.getPath(), "--bed", bed.getPath(), "--vcf", vcf.getPath(), "--output", output.getPath(), "--father", "Father", "--mother", "Mother", "-Z"};
      final int ret = cli.mainInit(args, out.outputStream(), err.printStream());
      assertEquals(err.toString(), 0, ret);
      final String result = FileUtils.fileToString(output);
      mNano.check("crossovers_annotated_repaired.vcf", result, true);
    }
  }
}
