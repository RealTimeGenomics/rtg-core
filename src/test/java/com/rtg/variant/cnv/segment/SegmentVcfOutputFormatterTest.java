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

import java.io.IOException;

import com.rtg.launcher.AbstractNanoTest;
import com.rtg.reader.ArraySequencesReader;
import com.rtg.util.StringUtils;
import com.rtg.util.TestUtils;

/**
 * Test the corresponding class.
 */
public class SegmentVcfOutputFormatterTest extends AbstractNanoTest {

  public void test() throws IOException {
    final SegmentVcfOutputFormatter formatter = new SegmentVcfOutputFormatter(new ArraySequencesReader(StringUtils.repeat("A", 1000), StringUtils.repeat("G", 1000), StringUtils.repeat("AG", 500)), 0.0, 1, "SAMPLE");
    final String header = formatter.header().toString();
    TestUtils.containsAll(header,
      "##contig=<ID=sequence 0,length=1000>",
      "##contig=<ID=sequence 1,length=1000>",
      "##contig=<ID=sequence 2,length=1000>",
      "##ALT=<ID=DEL,Description=\"Deletion\">",
      "##ALT=<ID=DUP,Description=\"Duplication\">",
      "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">",
      "##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise structural variation\">",
      "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">",
      "##INFO=<ID=BC,Number=1,Type=Integer,Description=\"Number of bins contained within the region\">",
      "##INFO=<ID=CIPOS,Number=2,Type=Integer,Description=\"Confidence interval around POS for imprecise variants\">",
      "##INFO=<ID=CIEND,Number=2,Type=Integer,Description=\"Confidence interval around END for imprecise variants\">",
      "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
      "#FORMAT=<ID=SQS,Number=1,Type=Float,Description=\"Segment quality score\">",
      "##FORMAT=<ID=RDR,Number=1,Type=Float,Description=\"Mean normalized RD ratio with respect to control\">",
      "##FORMAT=<ID=LR,Number=1,Type=Float,Description=\"Log2 of RD ratio with respect to control\">",
      //"##FORMAT=<ID=NSC,Number=1,Type=Float,Description=\"Mean normalized case coverage of the segment\">",
      //"##FORMAT=<ID=NCC,Number=1,Type=Float,Description=\"Mean normalized control coverage of the segment\">",
      "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE"
      );
    final String sb = String.valueOf(formatter.vcfRecord(null, new Segment("sequence 1", 1, 42, 4.3, 0, 1, 1), new Segment("sequence 1", 60, 70, 8.0, 48.0, 1, 1)))
      + '\n'
      + formatter.vcfRecord(new Segment("sequence 1", 1, 42, 4.3, 0, 1, 1), new Segment("sequence 1", 60, 70, 8.0, 48.0, 1, 1), null)
      + '\n';
    mNano.check("svof-example.txt", sb);
  }
}
