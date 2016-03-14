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

package com.rtg.variant.bayes.complex;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;

import com.rtg.reader.FastaUtils;
import com.rtg.util.MathUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.variant.NoQualityException;
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
    for (int i = 0; i < length; i++) {
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

  public void testHypothesisCorrection() throws IOException {
    final VariantOutputOptions params = VariantParams.builder().create();
    final VariantParams vp = HypothesesComplexTest.getVariantParams(0.5, 0.5, 0.1);
    final ComplexTemplate comtem = new ComplexTemplate(new byte[]{0, 0}, "foo", 0, 2);
    final HypothesesComplex hyp = HypothesesComplex.makeComplexHypotheses(comtem, Arrays.asList(match("AA", 5)), LogPossibility.SINGLETON, true, vp, null);
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
    for (int i = 0; i < sam.getReadString().length(); i++) {
      sb.append('`');
    }
    final int start = fixedLeft ? 5 : 0;
    return new AlignmentMatch(new VariantAlignmentRecord(sam), null, sam.getReadString(), FastaUtils.asciiToRawQuality(sb.toString()), 0, start, length, mapq, fixedLeft, fixedRight);
  }


  public void testIncompleteOverlap() throws IOException {
    final VariantParams vp = HypothesesComplexTest.getVariantParams(0.5, 0.5, 0.1);
    final ComplexTemplate comtem = new ComplexTemplate(new byte[]{0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, "foo", 0, 10);
    final HypothesesComplex hyp = HypothesesComplex.makeComplexHypotheses(comtem, Arrays.asList(partialMatch("AAAAAAAAAA", "10=", 5, true, true)), LogPossibility.SINGLETON, true, vp, null);
    final StatisticsComplex stat = new StatisticsComplex(hyp.description(), comtem.getLength());
    for (int i = 0; i < 10; i++) {
      stat.increment(new EvidenceComplex(hyp, partialMatch("AA", "2=", 5, true, false), comtem, vp, EvidenceComplexTest.getChooser()), 0);
    }
    assertEquals(2.000, stat.coverage(), 1e-3);
    assertEquals(0.632, stat.totalError(), 1e-3);
  }

  public void testIncompleteOverlapOtherEnd() throws IOException {
    final VariantParams vp = HypothesesComplexTest.getVariantParams(0.5, 0.5, 0.1);
    final ComplexTemplate comtem = new ComplexTemplate(new byte[]{0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, "foo", 0, 10);
    final HypothesesComplex hyp = HypothesesComplex.makeComplexHypotheses(comtem, Arrays.asList(partialMatch("AAAAAAAAAA", "10=", 5, true, true)), LogPossibility.SINGLETON, true, vp, null);
    final StatisticsComplex stat = new StatisticsComplex(hyp.description(), comtem.getLength());
    for (int i = 0; i < 10; i++) {
      stat.increment(new EvidenceComplex(hyp, partialMatch("AA", "2=", 5, false, true), comtem, vp, EvidenceComplexTest.getChooser()), 0);
    }
    assertEquals(2.000, stat.coverage(), 1e-3);
    assertEquals(0.632, stat.totalError(), 1e-3);
  }

  public void testExactCoverage() throws IOException {
    final VariantParams vp = HypothesesComplexTest.getVariantParams(0.5, 0.5, 0.1);
    final ComplexTemplate comtem = new ComplexTemplate(new byte[]{0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, "foo", 0, 10);
    final HypothesesComplex hyp = HypothesesComplex.makeComplexHypotheses(comtem, Arrays.asList(partialMatch("AAAAAAAAAA", "10=", 5, true, true)), LogPossibility.SINGLETON, true, vp, null);
    final StatisticsComplex stat = new StatisticsComplex(hyp.description(), comtem.getLength());
    for (int i = 0; i < 12; i++) {
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



  public void testCoverage() throws Exception {
    Diagnostic.setLogStream();
    final ArrayList<AlignmentMatch> ml = new ArrayList<>();
    for (int i = 0; i < 3; i++) {
      ml.add(match("AAA", 0));
    }
    for (int i = 0; i < 10; i++) {
      ml.add(match("GGG", 20));
    }
    final VariantParams vp = HypothesesComplexTest.getVariantParams(0.5, 0.5, 0.1);
    final ComplexTemplate cot = new ComplexTemplate(new byte[] {1, 2, 3}, "Foo", 0, 3);
    final HypothesesComplex hyp = HypothesesComplex.makeComplexHypotheses(cot, ml, LogPossibility.SINGLETON, true, vp, null);
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
    public double baseError(int index) throws NoQualityException, IndexOutOfBoundsException {
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
