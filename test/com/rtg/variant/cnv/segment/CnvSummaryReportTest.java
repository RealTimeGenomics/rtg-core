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

import org.junit.After;
import org.junit.Before;
import org.junit.Test;

import com.rtg.util.StringUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.intervals.RangeList;
import com.rtg.util.intervals.SimpleRangeMeta;
import com.rtg.util.intervals.ReferenceRanges;
import com.rtg.util.io.TestDirectory;
import com.rtg.util.test.FileHelper;
import com.rtg.util.test.NanoRegression;

/**
 */
public class CnvSummaryReportTest {
  protected NanoRegression mNano = null;

  @Before
  public void setUp() {
    Diagnostic.setLogStream();
    mNano = new NanoRegression(this.getClass(), false);
  }

  @After
  public void tearDown() throws Exception {
    try {
      mNano.finish();
    } finally {
      mNano = null;
    }
  }

  @Test
  public void test() throws IOException {
    try (TestDirectory temp = new TestDirectory()) {
      final File vcfFile = FileHelper.resourceToFile("com/rtg/variant/cnv/segment/resources/singleRecord.vcf", new File(temp, "singleRecord.vcf"));
      final File output = new File(temp, "out");
      final ReferenceRanges<String> regions = new ReferenceRanges<>(true);
      final RangeList<String> range = new RangeList<>(new SimpleRangeMeta<>(4580, 4690, "one"));
      regions.put("19", range);
      final CnvSummaryReport cnvSummaryReport = new CnvSummaryReport(regions);
      cnvSummaryReport.report(vcfFile, output);
      mNano.check("singleRecordOutput.txt", StringUtils.grepMinusV(FileHelper.fileToString(output), "^#[^c]"));
    }
  }


}
