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
package com.rtg.variant.util;

import java.util.Arrays;

import com.rtg.util.MathUtils;
import com.rtg.variant.VariantAlignmentRecord;
import com.rtg.variant.VariantParams;
import com.rtg.variant.VariantParamsBuilder;
import com.rtg.variant.util.arithmetic.LogApproximatePossibility;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMRecord;

import junit.framework.TestCase;

/**
 */
public class VariantUtilsTest extends TestCase {

  /**
   * Test method for {@link VariantUtils#phredToProb(int)}.
   */
  public final void testScoreToProb() {
    assertEquals(1.0, VariantUtils.phredToProb(0), 0.000001);
    assertEquals(0.316228, VariantUtils.phredToProb(5), 0.000001);
    assertEquals(0.1, VariantUtils.phredToProb(10), 0.000001);
    assertEquals(0.01, VariantUtils.phredToProb(20), 0.000001);
    assertEquals(0.0, VariantUtils.phredToProb(63), 0.000001);
    assertEquals(1e-50, VariantUtils.phredToProb(500), 1e-56);
}

  /**
   * Test method for {@link VariantUtils#phredToProb(int)}.
   */
  public final void testScoreToProbBad() {
    try {
      VariantUtils.phredToProb(-1);
      fail();
    } catch (final RuntimeException e) {
      assertTrue(e.getMessage().contains("-1"));
    }
//    try {
//      VariantUtils.phredToProb(255);
//      fail();
//    } catch (final RuntimeException e) {
//      assertTrue(e.getMessage().contains("255"));
//    }
  }

  public void testLogSum() {
    checkLogSum(1.0, 1.0);
    checkLogSum(0.0, 1.0);
    checkLogSum(0.0, 0.0);
    checkLogSum(0.01, 0.0000001);
    assertEquals(Double.POSITIVE_INFINITY, VariantUtils.logSum(Double.POSITIVE_INFINITY, 1.0), 0.000001);
    assertEquals(Double.POSITIVE_INFINITY, VariantUtils.logSum(1.0, Double.POSITIVE_INFINITY), 0.000001);
    assertEquals(Double.POSITIVE_INFINITY, VariantUtils.logSum(Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY), 0.000001);
  }

  private void checkLogSum(final double x, final double y) {
    assertEquals(Math.log(x + y), VariantUtils.logSum(Math.log(x), Math.log(y)), 0.000001);
    assertEquals(Math.log(x + y), VariantUtils.logSum(Math.log(y), Math.log(x)), 0.000001);
  }

  public void testLogSubtract() {
    checkLogSubtract(1.0, 1.0);
    assertEquals(0, Double.compare(Double.NaN, VariantUtils.logSubtract(/*Math.log(0.0)*/Double.NEGATIVE_INFINITY, /*Math.log(1.0)*/0.0)));
    checkLogSubtract(0.0, 0.0);
    checkLogSubtract(0.01, 0.0000001);
    assertEquals(Double.POSITIVE_INFINITY, VariantUtils.logSubtract(Double.POSITIVE_INFINITY, 1.0), 0.000001);
    assertEquals(0, Double.compare(Double.NaN, VariantUtils.logSubtract(1.0, Double.POSITIVE_INFINITY)));
    //assertEquals(Double.NaN, VariantUtils.logSubtract(Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY), 0.000001);
    assertEquals(0, Double.compare(Double.NaN, VariantUtils.logSubtract(Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY)));
  }

  private void checkLogSubtract(final double x, final double y) {
    assertEquals(Math.log(x - y), VariantUtils.logSubtract(Math.log(x), Math.log(y)), 0.000001);
  }


  public void testLogSumApproximation() {
    checkLogSumApproximation(1.0, 1.0);
    checkLogSumApproximation(0.0, 1.0);
    checkLogSumApproximation(0.0, 0.0);
    checkLogSumApproximation(0.01, 0.0000001);
    assertEquals(Double.POSITIVE_INFINITY, VariantUtils.logSumApproximation(Double.POSITIVE_INFINITY, 1.0), 0.000001);
    assertEquals(Double.POSITIVE_INFINITY, VariantUtils.logSumApproximation(1.0, Double.POSITIVE_INFINITY), 0.000001);
    assertEquals(Double.POSITIVE_INFINITY, VariantUtils.logSumApproximation(Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY), 0.000001);
  }

  private void checkLogSumApproximation(final double x, final double y) {
    assertEquals(Math.log(x + y), VariantUtils.logSumApproximation(Math.log(x), Math.log(y)), VariantUtils.MIN_X);
    assertEquals(Math.log(x + y), VariantUtils.logSumApproximation(Math.log(y), Math.log(x)), VariantUtils.MIN_X);
  }

  public void testIndexConversion() {
    for (int i = 0; i < VariantUtils.SIZE; i++) {
      checkIndexConversion(i + VariantUtils.START_INDEX);
    }
    assertEquals(VariantUtils.START_INDEX + VariantUtils.SIZE - 1, VariantUtils.doubleToIndex(VariantUtils.MAX_X));
  }

  private void checkIndexConversion(final int index) {
    final double value = VariantUtils.indexToDouble(index);
    final int index2 = VariantUtils.doubleToIndex(value);
    final double value2 = VariantUtils.indexToDouble(index2);
    assertEquals(index, index2);
    assertEquals(value, value2);
  }

  public void testLog1pExpApproximation() {
    assertEquals(VariantUtils.log1pExpApproximation(0.0), VariantUtils.log1pExpApproximation(VariantUtils.MIN_X - VariantUtils.MIN_X / 10.0));
    assertEquals(0.0, VariantUtils.log1pExpApproximation(VariantUtils.MAX_X + 1));
    assertEquals(VariantUtils.VALUES[0], VariantUtils.log1pExpApproximation(VariantUtils.MIN_X));
    assertEquals(VariantUtils.VALUES[3], VariantUtils.log1pExpApproximation(VariantUtils.indexToDouble(VariantUtils.START_INDEX + 3)));
    assertEquals(VariantUtils.VALUES[0], VariantUtils.logExpFunction((VariantUtils.indexToDouble(VariantUtils.START_INDEX) + VariantUtils.indexToDouble(VariantUtils.START_INDEX + 1)) / 2.0));
  }

  public void testLogExpFunction() {
    for (double i = 0; i < 20; i += 0.1) {
      assertEquals(Math.log1p(Math.exp(-i)), VariantUtils.logExpFunction(i));
    }
  }

  //simple case - get it from record
  public void testReadScoreFromSamRecord1() {
    final boolean mated = true;
    final VariantParamsBuilder vpb = new VariantParamsBuilder();
    final VariantParams params = vpb.create();
    final SAMRecord sam = new SAMRecord(new SAMFileHeader());
    sam.setReadPairedFlag(mated);
    sam.setProperPairFlag(mated);
    sam.setAttribute("NH", 1);
    sam.setMappingQuality(42);
    assertEquals(42, VariantUtils.readScoreFromAlignmentRecord(new VariantAlignmentRecord(sam), params));
  }

  //use max mated score to reduce
  public void testReadScoreFromSamRecord2() {
    final boolean mated = true;
    final VariantParamsBuilder vpb = new VariantParamsBuilder();
    vpb.maxMatedReadQuality(21);
    final VariantParams params = vpb.create();
    final SAMRecord sam = new SAMRecord(new SAMFileHeader());
    sam.setReadPairedFlag(mated);
    sam.setProperPairFlag(mated);
    sam.setAttribute("NH", 1);
    sam.setMappingQuality(42);
    assertEquals(21, VariantUtils.readScoreFromAlignmentRecord(new VariantAlignmentRecord(sam), params));
  }

  //Check that unmated doesnt affect simple case
  public void testReadScoreFromSamRecord3() {
    final boolean mated = false;
    final VariantParamsBuilder vpb = new VariantParamsBuilder();
    final VariantParams params = vpb.create();
    final SAMRecord sam = new SAMRecord(new SAMFileHeader());
    sam.setReadPairedFlag(mated);
    sam.setProperPairFlag(mated);
    sam.setAttribute("NH", 1);
    sam.setMappingQuality(42);
    assertEquals(42, VariantUtils.readScoreFromAlignmentRecord(new VariantAlignmentRecord(sam), params));
  }

  //unmated max clips score
  public void testReadScoreFromSamRecord4() {
    final boolean mated = false;
    final VariantParamsBuilder vpb = new VariantParamsBuilder();
    vpb.maxUnmatedReadQuality(21);
    final VariantParams params = vpb.create();
    final SAMRecord sam = new SAMRecord(new SAMFileHeader());
    sam.setReadPairedFlag(mated);
    sam.setProperPairFlag(mated);
    sam.setAttribute("NH", 1);
    sam.setMappingQuality(42);
    assertEquals(21, VariantUtils.readScoreFromAlignmentRecord(new VariantAlignmentRecord(sam), params));
  }

  //NH=1 value forces result - use mated default
  public void testReadScoreFromSamRecord5() {
    final boolean mated = true;
    final VariantParamsBuilder vpb = new VariantParamsBuilder();
    vpb.defaultMatedReadQuality(23);
    final VariantParams params = vpb.create();
    final SAMRecord sam = new SAMRecord(new SAMFileHeader());
    sam.setReadPairedFlag(mated);
    sam.setProperPairFlag(mated);
    sam.setAttribute("NH", 1);
    sam.setMappingQuality(255);
    assertEquals(23, VariantUtils.readScoreFromAlignmentRecord(new VariantAlignmentRecord(sam), params));
  }

  //NH=1 value forces result - use unmated default
  public void testReadScoreFromSamRecord6() {
    final boolean mated = false;
    final VariantParamsBuilder vpb = new VariantParamsBuilder();
    vpb.defaultUnmatedReadQuality(25);
    final VariantParams params = vpb.create();
    final SAMRecord sam = new SAMRecord(new SAMFileHeader());
    sam.setReadPairedFlag(mated);
    sam.setProperPairFlag(mated);
    sam.setAttribute("NH", 1);
    sam.setMappingQuality(255);
    assertEquals(25, VariantUtils.readScoreFromAlignmentRecord(new VariantAlignmentRecord(sam), params));
  }

  //NH=null value forces result - use mated default
  public void testReadScoreFromSamRecord7() {
    final boolean mated = true;
    final VariantParamsBuilder vpb = new VariantParamsBuilder();
    vpb.defaultMatedReadQuality(23);
    final VariantParams params = vpb.create();
    final SAMRecord sam = new SAMRecord(new SAMFileHeader());
    sam.setReadPairedFlag(mated);
    sam.setProperPairFlag(mated);
    sam.setAttribute("NH", 1);
    sam.setMappingQuality(255);
    assertEquals(23, VariantUtils.readScoreFromAlignmentRecord(new VariantAlignmentRecord(sam), params));
  }

  //NH=1 value forces result - use mated default
  public void testReadScoreFromSamRecord8() {
    final boolean mated = true;
    final VariantParamsBuilder vpb = new VariantParamsBuilder();
    vpb.defaultMatedReadQuality(23);
    final VariantParams params = vpb.create();
    final SAMRecord sam = new SAMRecord(new SAMFileHeader());
    sam.setReadPairedFlag(mated);
    sam.setProperPairFlag(mated);
    sam.setMappingQuality(255);
    assertEquals(23, VariantUtils.readScoreFromAlignmentRecord(new VariantAlignmentRecord(sam), params));
  }

  //NH=null value forces result - use unmated default
  public void testReadScoreFromSamRecord9() {
    final boolean mated = false;
    final VariantParamsBuilder vpb = new VariantParamsBuilder();
    vpb.defaultUnmatedReadQuality(25);
    final VariantParams params = vpb.create();
    final SAMRecord sam = new SAMRecord(new SAMFileHeader());
    sam.setReadPairedFlag(mated);
    sam.setProperPairFlag(mated);
    sam.setMappingQuality(255);
    assertEquals(25, VariantUtils.readScoreFromAlignmentRecord(new VariantAlignmentRecord(sam), params));
  }

  //NH=2
  public void testReadScoreFromSamRecord10() {
    final VariantParamsBuilder vpb = new VariantParamsBuilder();
    final VariantParams params = vpb.create();
    final SAMRecord sam = new SAMRecord(new SAMFileHeader());
    sam.setAttribute("NH", 2);
    sam.setMappingQuality(255);
    assertEquals(VariantUtils.phredFromN(2), VariantUtils.readScoreFromAlignmentRecord(new VariantAlignmentRecord(sam), params));
  }

  //NH=0 ?how can this ever happen
  public void testReadScoreFromSamRecord11() {
    final VariantParamsBuilder vpb = new VariantParamsBuilder();
    final VariantParams params = vpb.create();
    final SAMRecord sam = new SAMRecord(new SAMFileHeader());
    sam.setAttribute("NH", 0);
    sam.setMappingQuality(255);
    assertEquals(0, VariantUtils.readScoreFromAlignmentRecord(new VariantAlignmentRecord(sam), params));
  }

  //ignore actual quality use mated default
  public void testReadScoreFromSamRecord12() {
    final boolean mated = true;
    final VariantParamsBuilder vpb = new VariantParamsBuilder();
    vpb.defaultMatedReadQuality(17);
    vpb.ignoreReadQuality(true);
    final VariantParams params = vpb.create();
    final SAMRecord sam = new SAMRecord(new SAMFileHeader());
    sam.setReadPairedFlag(mated);
    sam.setProperPairFlag(mated);
    sam.setAttribute("NH", 1);
    sam.setMappingQuality(42);
    assertEquals(17, VariantUtils.readScoreFromAlignmentRecord(new VariantAlignmentRecord(sam), params));
  }

  //ignore actual quality use unmated default
  public void testReadScoreFromSamRecord13() {
    final boolean mated = false;
    final VariantParamsBuilder vpb = new VariantParamsBuilder();
    vpb.defaultUnmatedReadQuality(19);
    vpb.ignoreReadQuality(true);
    final VariantParams params = vpb.create();
    final SAMRecord sam = new SAMRecord(new SAMFileHeader());
    sam.setReadPairedFlag(mated);
    sam.setProperPairFlag(mated);
    sam.setAttribute("NH", 1);
    sam.setMappingQuality(42);
    assertEquals(19, VariantUtils.readScoreFromAlignmentRecord(new VariantAlignmentRecord(sam), params));
  }

  //ignore NH=2
  public void testReadScoreFromSamRecord14() {
    final boolean mated = true;
    final VariantParamsBuilder vpb = new VariantParamsBuilder();
    vpb.defaultMatedReadQuality(17);
    vpb.ignoreReadQuality(true);
    final VariantParams params = vpb.create();
    final SAMRecord sam = new SAMRecord(new SAMFileHeader());
    sam.setReadPairedFlag(mated);
    sam.setProperPairFlag(mated);
    sam.setAttribute("NH", 2);
    sam.setMappingQuality(42);
    assertEquals(VariantUtils.phredFromN(2), VariantUtils.readScoreFromAlignmentRecord(new VariantAlignmentRecord(sam), params));
  }

  private static final double PHRED_SCALE = -10.0 / Math.log(10.0);

  private static int phredCheck(final int n) {
    final double x = -1.0 / n;
    final long phred = (int) MathUtils.round(PHRED_SCALE * Math.log1p(x));
    if (phred > 63) {
      return 63;
    }
    return (int) phred;
  }

  private void checkPhred(final int phred, final int n) {
    assertEquals(phred, VariantUtils.phredFromN(n));
    assertEquals(phredCheck(n), VariantUtils.phredFromN(n));
  }

  public void testPhred() {
    checkPhred(3, 2);
    checkPhred(2, 3);
    checkPhred(1, 4);
    checkPhred(1, 5);
    checkPhred(1, 6);
    checkPhred(1, 7);
    checkPhred(1, 8);
    checkPhred(1, 9);
    checkPhred(0, 10);
    try {
      VariantUtils.phredFromN(1);
    } catch (final RuntimeException e) {
      assertEquals("Invalid record count:1", e.getMessage());
    }
  }

  public void testNormalisation() {
    final double[] normed = VariantUtils.normalisePossibilities(new double[]{0.4565873371437609, -20.351968280794335, -36.661425482130056, -36.62905520787501, -37.43049091089874, -19.82279123723306, -53.125023572451596, -69.40211049953227, -70.17117592830095, -36.13224843856878, -53.09265329819655, -70.203546202556, -36.099878164313736, -53.89408900122028, -36.90131386733746}, LogApproximatePossibility.SINGLETON);
    assertTrue(Arrays.equals(new double[] {0.0, -20.808555617938097, -37.118012819273815, -37.08564254501877, -37.8870782480425, -20.27937857437682, -53.581610909595355, -69.85869783667603, -70.62776326544471, -36.58883577571254, -53.54924063534031, -70.66013353969976, -36.556465501457495, -54.35067633836404, -37.357901204481216}, normed));
  }

  public void testNormalizePair() {
    checkNormalizePair("", "", "");
    checkNormalizePair(":A", "", "A");
    checkNormalizePair("C:A", "C", "A");
    checkNormalizePair("A", "A", "A");
    checkNormalizePair("A:C", "A", "C");
    try {
      checkNormalizePair("C:A:A", "C", "A");
      fail();
    } catch (final IllegalArgumentException e) {
      assertEquals("Invalid variation: C:A:A", e.getMessage());
    }
  }

  private void checkNormalizePair(final String variation, final String left, final String right) {
    final String[] res = VariantUtils.normalizePair(variation);
    assertEquals(2, res.length);
    assertEquals(left, res[0]);
    assertEquals(right, res[1]);
  }
}
