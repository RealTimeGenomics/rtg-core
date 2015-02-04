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

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMRecord;

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
    for (int i = 0; i < length; i++) {
      sb.append('`');
    }
    return new AlignmentMatch(new VariantAlignmentRecord(sam), null, ins, FastaUtils.asciiToRawQuality(sb.toString()), 0, 0, length, mapq, fixedLeft, fixedRight);
  }
  static AlignmentMatch match(final String ins, int mapq) {
    return match(ins, mapq, true, true);
  }
  public void testTotals() throws Exception {
    Diagnostic.setLogStream();
    final ArrayList<AlignmentMatch> ml = new ArrayList<>();
    for (int i = 0; i < 3; i++) {
      ml.add(match("A", 0));
    }
    for (int i = 0; i < 10; i++) {
      ml.add(match("G", 20));
    }
    final VariantParams vp = HypothesesComplexTest.getVariantParams(0.5, 0.5, 0.1);
    final ComplexTemplate cot = new ComplexTemplate(new byte[] {}, "", 0, 0);
    final ArrayList<AlignmentMatch> hypMatches = new ArrayList<>();
    hypMatches.add(match("T", 7));
    hypMatches.addAll(ml);
    final HypothesesComplex hyp = HypothesesComplex.makeComplexHypotheses(cot, hypMatches, LogPossibility.SINGLETON, true, vp, null);
    final StatisticsComplex cmpx = new StatisticsComplex(hyp.description(), cot.getLength());
    final ModelInterface<?> m = new Model<>(hyp, cmpx);
    for (final AlignmentMatch match : ml) {
      m.increment(new EvidenceComplex(hyp, match, cot, vp, EvidenceComplexTest.getChooser()));
    }
    final VariantOutputOptions params = VariantParams.builder().create();
    final VariantSample vs = new VariantSample(null, null, true, null, VariantSample.DeNovoStatus.UNSPECIFIED, null);
    cmpx.addCountsToSample(vs, null, params);
    assertEquals(3 / 13.0, cmpx.ambiguityRatio());
  }
}
