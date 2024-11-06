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
