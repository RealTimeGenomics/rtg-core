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
package com.rtg.variant.bayes;

import java.util.ArrayList;

import com.rtg.reader.FastaUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.variant.VariantAlignmentRecord;
import com.rtg.variant.VariantOutputOptions;
import com.rtg.variant.VariantParams;
import com.rtg.variant.VariantSample;
import com.rtg.variant.bayes.complex.ComplexTemplate;
import com.rtg.variant.bayes.complex.EvidenceComplex;
import com.rtg.variant.bayes.complex.EvidenceComplexTest;
import com.rtg.variant.bayes.complex.HypothesesComplex;
import com.rtg.variant.bayes.complex.HypothesesComplexTest;
import com.rtg.variant.bayes.complex.StatisticsComplex;
import com.rtg.variant.match.AlignmentMatch;
import com.rtg.variant.util.arithmetic.LogPossibility;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;

/**
 */
public class StatisticsDoubleTest extends AbstractStatisticsTest {

  @Override
  protected Statistics<?> getStatistics(Description d) {
    return new StatisticsDouble(d);
  }

  static AlignmentMatch match(final String ins, int mapq, boolean fixedLeft, boolean fixedRight) {
    final int length = ins.length();
    final SAMRecord sam = new SAMRecord(new SAMFileHeader());
    sam.setAlignmentStart(1);
    sam.setReadString("CCCCC" + ins + "GGGGG");
    sam.setCigarString("5=" + length + "I5=");
    final StringBuilder sb = new StringBuilder();
    for (int i = 0; i < length; ++i) {
      sb.append('`');
    }
    return new AlignmentMatch(new VariantAlignmentRecord(sam), null, ins, FastaUtils.asciiToRawQuality(sb.toString()), 0, 0, length, mapq, fixedLeft, fixedRight);
  }

  static AlignmentMatch match(final String ins, int mapq) {
    return match(ins, mapq, true, true);
  }

  public void testTotals() {
    Diagnostic.setLogStream();
    final ArrayList<AlignmentMatch> ml = new ArrayList<>();
    for (int i = 0; i < 3; ++i) {
      ml.add(match("A", 0));
    }
    for (int i = 0; i < 10; ++i) {
      ml.add(match("G", 20));
    }
    final VariantParams vp = HypothesesComplexTest.getVariantParams(0.5, 0.5, 0.1);
    final ComplexTemplate cot = new ComplexTemplate(new byte[1], "A", 0, 1);
    final ArrayList<AlignmentMatch> hypMatches = new ArrayList<>();
    hypMatches.add(match("T", 7));
    hypMatches.addAll(ml);
    cot.setComplexContext(HypothesesComplex.createComplexDescription(hypMatches, cot, null, vp.pruneHypotheses(), vp.maxComplexHypotheses()), LogPossibility.SINGLETON);
    final HypothesesComplex hyp = HypothesesComplex.makeComplexHypotheses(cot, true, vp);
    final StatisticsComplex cmpx = new StatisticsComplex(hyp.description(), cot.getLength());
    final ModelInterface<?> m = new Model<>(hyp, cmpx, new NoAlleleBalance());
    for (final AlignmentMatch match : ml) {
      m.increment(new EvidenceComplex(hyp, match, cot, vp, EvidenceComplexTest.getChooser()));
    }
    final VariantOutputOptions params = VariantParams.builder().create();
    final VariantSample vs = new VariantSample(null, null, true, null, VariantSample.DeNovoStatus.UNSPECIFIED, null);
    cmpx.addCountsToSample(vs, null, params);
    assertEquals(3 / 13.0, cmpx.ambiguityRatio());
  }
}
