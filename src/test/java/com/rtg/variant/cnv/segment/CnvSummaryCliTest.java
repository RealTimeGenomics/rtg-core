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

package com.rtg.variant.cnv.segment;

import java.io.File;
import java.io.IOException;

import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.AbstractCliTest;
import com.rtg.launcher.MainResult;
import com.rtg.util.StringUtils;
import com.rtg.util.TestUtils;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.TestDirectory;
import com.rtg.util.test.FileHelper;

/**
 * Tests the corresponding class.
 */
public class CnvSummaryCliTest extends AbstractCliTest {

  @Override
  protected AbstractCli getCli() {
    return new CnvSummaryCli();
  }

  public void testHelp() {
    checkHelp("intersection of CNVs",
      "VCF file containing CNV",
      "BED output file",
      "BED file supplying gene-scale",
      "also report no alteration regions",
      "Reporting",
      "File Input/Output",
      "do not gzip the output"
      );
  }

  public void testErrors() {
    TestUtils.containsAll(checkHandleFlagsErr(), "Error: You must provide values for -i FILE -o FILE --summary-regions");
  }

  public void testActualRun() throws IOException {
    try (final TestDirectory dir = new TestDirectory()) {
      final File vcf = FileHelper.resourceToFile("com/rtg/variant/cnv/segment/resources/my.vcf", new File(dir, "my.vcf"));
      final File regions = FileHelper.resourceToFile("com/rtg/variant/cnv/segment/resources/regions.bed", new File(dir, "regions.bed"));
      final File output = new File(dir, "output.bed");
      final MainResult result = MainResult.run(getCli(), "-i", vcf.getPath(), "-o", output.getPath(), "--summary-regions", regions.getPath(), "-Z");
      assertEquals(0, result.rc());
      mNano.check("expected.summary.bed", StringUtils.grepMinusV(FileUtils.fileToString(output), "^#"));
      mNano.check("expected.sum-out.txt", result.out());
    }

  }

}
