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
package com.rtg.variant.avr;

import java.io.File;

import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.AbstractCliTest;
import com.rtg.util.TestUtils;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.util.io.TestDirectory;
import com.rtg.util.test.FileHelper;

/**
 */
public class PredictCliTest extends AbstractCliTest {

  @Override
  protected AbstractCli getCli() {
    return new PredictCli();
  }

  public void test() throws Exception {
    try (final TestDirectory dir = new TestDirectory()) {
      final File model = FileHelper.resourceToFile("com/rtg/variant/avr/resources/mlmodel.avr", new File(dir, "mlmodel.avr"));
      final File vcf = FileHelper.resourceToFile("com/rtg/variant/avr/resources/multisample.vcf", new File(dir, "multisample.vcf"));
      final File output = new File(dir, "output.vcf.gz");
      final MemoryPrintStream mps = new MemoryPrintStream();
      final int code = new PredictCli().mainInit(new String[] {"--avr-model", model.getPath(), "-i", vcf.getPath(), "-o", output.getPath()}, mps.outputStream(), mps.printStream());
      assertEquals(mps.toString(), 0, code);
      mNano.check("predict.vcf", TestUtils.sanitizeVcfHeader(FileHelper.gzFileToString(output)));
    }
  }

  public void testSample() throws Exception {
    try (final TestDirectory dir = new TestDirectory()) {
      final File model = FileHelper.resourceToFile("com/rtg/variant/avr/resources/mlmodel.avr", new File(dir, "mlmodel.avr"));
      final File vcf = FileHelper.resourceToFile("com/rtg/variant/avr/resources/multisample.vcf", new File(dir, "multisample.vcf"));
      final File output = new File(dir, "output.vcf.gz");
      final MemoryPrintStream mps = new MemoryPrintStream();
      final int code = new PredictCli().mainInit(new String[] {"-s", "NA12892", "--avr-model", model.getPath(), "-i", vcf.getPath(), "-o", output.getPath()}, mps.outputStream(), mps.printStream());
      assertEquals(mps.toString(), 0, code);
      mNano.check("predictSample.vcf", TestUtils.sanitizeVcfHeader(FileHelper.gzFileToString(output)));
    }
  }
}
