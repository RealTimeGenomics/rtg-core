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
      "do not gzip the output"
      );
  }

  public void testErrors() throws IOException {
    final String res = checkHandleFlagsErr();
    final String exp = getCFlags().getUsageHeader();
    assertTrue(res.contains(exp));
    TestUtils.containsAll(res, "Error: You must provide values for -i FILE -o FILE --summary-regions");
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
