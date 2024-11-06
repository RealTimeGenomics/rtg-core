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
