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

package com.rtg.variant.bayes.complex;

import java.util.ArrayList;
import java.util.Arrays;

import com.rtg.reader.FastaUtils;
import com.rtg.util.MathUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.variant.VariantAlignmentRecord;
import com.rtg.variant.VariantOutputOptions;
import com.rtg.variant.VariantParams;
import com.rtg.variant.VariantSample;
import com.rtg.variant.match.AlignmentMatch;
import com.rtg.variant.match.Match;
import com.rtg.variant.util.arithmetic.LogPossibility;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import junit.framework.TestCase;

/**
 */
public class StatisticsComplexTest extends TestCase {

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

  public void testEmpty() {
    final StatisticsComplex cmpx = new StatisticsComplex(new DescriptionComplex(Arrays.asList(match("AA", 5))), 2);
    final VariantOutputOptions params2 = VariantParams.builder().create();
    final VariantSample vs = new VariantSample(null, null, true, null, VariantSample.DeNovoStatus.UNSPECIFIED, null);
    cmpx.addCountsToSample(vs, null, params2);
    assertEquals(0.0, vs.getCorrection(), 1e-13);
    assertEquals(null, cmpx.ambiguityRatio());
  }

  public void testHypothesisCorrection() {
    final VariantOutputOptions params = VariantParams.builder().create();
    final VariantParams vp = HypothesesComplexTest.getVariantParams(0.5, 0.5, 0.1);
    final ComplexTemplate comtem = new ComplexTemplate(new byte[]{0, 0}, "foo", 0, 2);
    comtem.setComplexContext(HypothesesComplex.createComplexDescription(Arrays.asList(match("AA", 5)), comtem, null, vp.pruneHypotheses(), vp.maxComplexHypotheses()), LogPossibility.SINGLETON);
    final HypothesesComplex hyp = HypothesesComplex.makeComplexHypotheses(comtem, true, vp);
    final StatisticsComplex stat = new StatisticsComplex(hyp.description(), comtem.getLength());
    stat.increment(new EvidenceComplex(hyp, match("A", 20), comtem, vp, EvidenceComplexTest.getChooser()), 0);
    stat.increment(new EvidenceComplex(hyp, match("C", 21), comtem, vp, EvidenceComplexTest.getChooser()), 0);
    stat.increment(new EvidenceComplex(hyp, match("G", 22), comtem, vp, EvidenceComplexTest.getChooser()), 0);
    stat.increment(new EvidenceComplex(hyp, match("T", 23), comtem, vp, EvidenceComplexTest.getChooser()), 0);
    final VariantSample vs = new VariantSample(null, null, true, null, VariantSample.DeNovoStatus.UNSPECIFIED, null);
    stat.addCountsToSample(vs, null, params);
    //  A 1 0.010 C 1 0.008 G 1 0.006 T 1 0.005
    assertEquals(0.029, vs.getCorrection(), 1e-3);
  }

  static AlignmentMatch partialMatch(final String bases, String partialCigar, int mapq, boolean fixedLeft, boolean fixedRight) {
    final int length = bases.length();
    final SAMRecord sam = new SAMRecord(new SAMFileHeader());
    sam.setAlignmentStart(1);
    sam.setReadString((fixedLeft ? "CCCCC" : "") + bases + (fixedRight ? "GGGGG" : ""));
    sam.setCigarString((fixedLeft ? "5=" : "") + partialCigar + (fixedRight ? "5=" : ""));
    final StringBuilder sb = new StringBuilder();
    for (int i = 0; i < sam.getReadString().length(); ++i) {
      sb.append('`');
    }
    final int start = fixedLeft ? 5 : 0;
    return new AlignmentMatch(new VariantAlignmentRecord(sam), null, sam.getReadString(), FastaUtils.asciiToRawQuality(sb.toString()), 0, start, length, mapq, fixedLeft, fixedRight);
  }


  public void testIncompleteOverlap() {
    final VariantParams vp = HypothesesComplexTest.getVariantParams(0.5, 0.5, 0.1);
    final ComplexTemplate comtem = new ComplexTemplate(new byte[]{0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, "foo", 0, 10);
    comtem.setComplexContext(HypothesesComplex.createComplexDescription(Arrays.asList(partialMatch("AAAAAAAAAA", "10=", 5, true, true)), comtem, null, vp.pruneHypotheses(), vp.maxComplexHypotheses()), LogPossibility.SINGLETON);
    final HypothesesComplex hyp = HypothesesComplex.makeComplexHypotheses(comtem, true, vp);
    final StatisticsComplex stat = new StatisticsComplex(hyp.description(), comtem.getLength());
    for (int i = 0; i < 10; ++i) {
      stat.increment(new EvidenceComplex(hyp, partialMatch("AA", "2=", 5, true, false), comtem, vp, EvidenceComplexTest.getChooser()), 0);
    }
    assertEquals(2.000, stat.coverage(), 1e-3);
    assertEquals(0.632, stat.totalError(), 1e-3);
  }

  public void testIncompleteOverlapOtherEnd() {
    final VariantParams vp = HypothesesComplexTest.getVariantParams(0.5, 0.5, 0.1);
    final ComplexTemplate comtem = new ComplexTemplate(new byte[]{0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, "foo", 0, 10);
    comtem.setComplexContext(HypothesesComplex.createComplexDescription(Arrays.asList(partialMatch("AAAAAAAAAA", "10=", 5, true, true)), comtem, null, vp.pruneHypotheses(), vp.maxComplexHypotheses()), LogPossibility.SINGLETON);
    final HypothesesComplex hyp = HypothesesComplex.makeComplexHypotheses(comtem, true, vp);
    final StatisticsComplex stat = new StatisticsComplex(hyp.description(), comtem.getLength());
    for (int i = 0; i < 10; ++i) {
      stat.increment(new EvidenceComplex(hyp, partialMatch("AA", "2=", 5, false, true), comtem, vp, EvidenceComplexTest.getChooser()), 0);
    }
    assertEquals(2.000, stat.coverage(), 1e-3);
    assertEquals(0.632, stat.totalError(), 1e-3);
  }

  public void testExactCoverage() {
    final VariantParams vp = HypothesesComplexTest.getVariantParams(0.5, 0.5, 0.1);
    final ComplexTemplate comtem = new ComplexTemplate(new byte[]{0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, "foo", 0, 10);
    comtem.setComplexContext(HypothesesComplex.createComplexDescription(Arrays.asList(partialMatch("AAAAAAAAAA", "10=", 5, true, true)), comtem, null, vp.pruneHypotheses(), vp.maxComplexHypotheses()), LogPossibility.SINGLETON);
    final HypothesesComplex hyp = HypothesesComplex.makeComplexHypotheses(comtem, true, vp);
    final StatisticsComplex stat = new StatisticsComplex(hyp.description(), comtem.getLength());
    for (int i = 0; i < 12; ++i) {
      stat.increment(new EvidenceComplex(hyp, partialMatch("AA", "2=", 5, true, false), comtem, vp, EvidenceComplexTest.getChooser()), 0);
    }
    assertEquals(2.400, stat.exactCoverage(), 1e-3);
  }

  /*
  public void testAmbiguityCount() {
    final StatisticsComplex cmpx = new StatisticsComplex(new DescriptionComplex(Arrays.asList(match("AA", 5))), 0);
    final VariantOutputOptions params = VariantParams.builder().statistics(false).create();
    cmpx.incrementComplex(match("A", 0));
    cmpx.incrementComplex(match("G", 20));
    final VariantSample vs = new VariantSample(null, null, true, null, VariantSample.DeNovoStatus.UNSPECIFIED, null);
    cmpx.addCountsToSample(vs, null, params);
  }*/

  public void testStupidRounding() {
    //Java 1.7 has a special case in Math.round() for 0.49999999999999994
    //Dealing with by using our own copy of the Java 1.6 implementation.
    final String ref = "ACGTATGAATGCATGCATGCATGCTAGC"; // 28 long
    final StatisticsComplex complexStat = new StatisticsComplex(new DescriptionComplex(Arrays.asList(match(ref, 5))), 28);
    double cov = 0.0;
    cov += complexStat.coverageIncrement(new DummyMatch(4));
    cov += complexStat.coverageIncrement(new DummyMatch(6));
    cov += complexStat.coverageIncrement(new DummyMatch(4));
    assertEquals(1, MathUtils.round(cov));
  }



  public void testCoverage() {
    Diagnostic.setLogStream();
    final ArrayList<AlignmentMatch> ml = new ArrayList<>();
    for (int i = 0; i < 3; ++i) {
      ml.add(match("AAA", 0));
    }
    for (int i = 0; i < 10; ++i) {
      ml.add(match("GGG", 20));
    }
    final VariantParams vp = HypothesesComplexTest.getVariantParams(0.5, 0.5, 0.1);
    final ComplexTemplate cot = new ComplexTemplate(new byte[] {1, 2, 3}, "Foo", 0, 3);
    cot.setComplexContext(HypothesesComplex.createComplexDescription(ml, cot, null, vp.pruneHypotheses(), vp.maxComplexHypotheses()), LogPossibility.SINGLETON);
    final HypothesesComplex hyp = HypothesesComplex.makeComplexHypotheses(cot, true, vp);
    final StatisticsComplex cmpx = new StatisticsComplex(hyp.description(), cot.getLength());
    double cov = 0.0;
    for (final AlignmentMatch match : ml) {
      cov += cmpx.coverageIncrement(match);
    }
    assertEquals(13L, MathUtils.round(cov));
    cov += cmpx.coverageIncrement(match("AA", 20, true, false));
    assertEquals(14L, MathUtils.round(cov));
    cov += cmpx.coverageIncrement(match("A", 20, true, false));
    assertEquals(14L, MathUtils.round(cov));
    cov += cmpx.coverageIncrement(match("A", 20, true, false));
    assertEquals(14L, MathUtils.round(cov));
    cov += cmpx.coverageIncrement(match("A", 20, true, false));
    assertEquals(15L, MathUtils.round(cov));
  }

  private static final class DummyMatch extends Match {
    final int mLen;
    DummyMatch(int len) {
      mLen = len;
    }
    @Override
    public String readString() {
      return null;
    }
    @Override
    public int read(int index) {
      return 0;
    }
    @Override
    public double baseError(int index) {
      return 0;
    }
    @Override
    public double baseError() {
      return 0;
    }
    @Override
    public double mapError() {
      return 0;
    }
    @Override
    public int length() {
      return mLen;
    }
    @Override
    public boolean isFixedRight() {
      return false;
    }
    @Override
    public boolean isFixedLeft() {
      return false;
    }
  }

}
