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

import java.io.IOException;

import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.AbstractCliTest;
import com.rtg.util.TestUtils;

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
      "do not gzip the output"
      );
  }

  public void testErrors() throws IOException {
    final String res = checkHandleFlagsErr();
    final String exp = getCFlags().getUsageHeader();
    assertTrue(res.contains(exp));
    TestUtils.containsAll(res, "Error: You must provide values for -i FILE -o FILE --summary-regions");
  }

}
