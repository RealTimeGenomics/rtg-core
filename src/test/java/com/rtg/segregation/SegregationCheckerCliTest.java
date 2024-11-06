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

package com.rtg.segregation;

import java.io.File;
import java.io.IOException;

import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.AbstractCliTest;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.reference.ReferenceGenome;
import com.rtg.util.TestUtils;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.TestDirectory;
import com.rtg.util.test.FileHelper;

public class SegregationCheckerCliTest extends AbstractCliTest {

  @Override
  protected AbstractCli getCli() {
    return new SegregationCheckerCli();
  }

  public void testNano() throws IOException {
    try (final TestDirectory dir = new TestDirectory()) {
      final File template = ReaderTestUtils.getDNADir(">1\nACGT\n>2\nACGT\n>3\nACGT", new File(dir, "template"));
      FileUtils.stringToFile("version 1\neither def diploid linear", new File(template, ReferenceGenome.REFERENCE_FILE));
      final File vcf = FileHelper.resourceToFile("com/rtg/segregation/resources/crossover.vcf", new File(dir, "vcf.vcf"));
      final File bed = FileHelper.resourceToFile("com/rtg/segregation/resources/regions.bed", new File(dir, "regions.bed"));
      final File output = new File(dir, "out.vcf");
      checkMainInitOk("--template", template.getPath(), "--bed", bed.getPath(), "--vcf", vcf.getPath(), "--output", output.getPath(), "--father", "Father", "--mother", "Mother", "-Z");
      final String result = TestUtils.sanitizeVcfHeader(FileUtils.fileToString(output));
      mNano.check("crossovers_annotated.vcf", result, true);
    }
  }

  public void testNanoRepair() throws IOException {
    try (final TestDirectory dir = new TestDirectory()) {
      final File template = ReaderTestUtils.getDNADir(">1\nACGT\n>2\nACGT\n>3\nACGT", new File(dir, "template"));
      FileUtils.stringToFile("version 1\neither def diploid linear", new File(template, ReferenceGenome.REFERENCE_FILE));
      final File vcf = FileHelper.resourceToFile("com/rtg/segregation/resources/crossover.vcf", new File(dir, "vcf.vcf"));
      final File bed = FileHelper.resourceToFile("com/rtg/segregation/resources/regions.bed", new File(dir, "regions.bed"));
      final File output = new File(dir, "out.vcf");
      checkMainInitOk("--repair", "--template", template.getPath(), "--bed", bed.getPath(), "--vcf", vcf.getPath(), "--output", output.getPath(), "--father", "Father", "--mother", "Mother", "-Z");
      final String result = TestUtils.sanitizeVcfHeader(FileUtils.fileToString(output));
      mNano.check("crossovers_annotated_repaired.vcf", result, true);
    }
  }

}
