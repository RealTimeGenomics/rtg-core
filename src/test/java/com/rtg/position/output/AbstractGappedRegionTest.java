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
package com.rtg.position.output;


import static com.rtg.position.output.GapProbabilitiesScorer.PROB_FORMAT;

import java.io.IOException;

import com.rtg.launcher.BuildParams;
import com.rtg.mode.SequenceMode;
import com.rtg.util.array.WrappedIntArray;

import junit.framework.TestCase;

/**
 * Test GappedRegion
 * @param <T> enforce a type of gapped region
 */
public abstract class AbstractGappedRegionTest<T extends AbstractGappedRegion<T>> extends TestCase {

  protected final T getRegion(final int id, final int stepSize, final int wordSize) throws IOException {
    final GappedDistribution distr = new GappedDistribution(stepSize, wordSize, GappedDistribution.distrParams(3 * stepSize));
    final GapScorer prob = distr.probabilities();
    //System.err.println(prob);
    return getRegion(id, stepSize, wordSize, prob, new WrappedIntArray(new int[] {42, 42}));
  }

  protected final T getRegion(final int id, final int stepSize, final int wordSize, final GapScorer prob) throws IOException {
    //System.err.println(prob);
    return getRegion(id, stepSize, wordSize, prob, new WrappedIntArray(new int[] {42, 42}));
  }

  protected T getRegion(final int id, final int stepSize, final int wordSize, final GapScorer prob, final WrappedIntArray buildLengths) throws IOException {
//    final ReaderParams rp = new MockReaderParams(100, 1, SequenceMode.UNIDIRECTIONAL);
//    final SequenceParams sp = new SequenceParams(rp);
    final BuildParams params = BuildParams.builder().windowSize(wordSize).stepSize(stepSize).sequences(null).create();
    return getRegion(id, params, prob, buildLengths);
  }

  abstract T getRegion(final int id, final BuildParams params, final GapScorer prob, final WrappedIntArray buildLengths);
  GapBuckets<T> makeBuckets(final int stepSize) throws IOException {
    final BuildParams buildParams = GapBucketsTest.makeParams(stepSize, new int[] {30, 31, 48}, 0, SequenceMode.UNIDIRECTIONAL);
    final GapBucketsInfo bucketInfo = new GapBucketsInfo(buildParams, 1, 16, 1000);
    return new GapBuckets<>(bucketInfo);
  }

  protected void check(final String expected1, final String expected2) throws IOException {
    final GapBuckets<T> bu = makeBuckets(8);
    final T gr = getRegion(1, 8, 12);
    final T gr2 = getRegion(2, 8, 12);

    gr.integrity();
    assertFalse(gr.isInitialized());
    assertEquals("empty", gr.toString());
    assertEquals(1, gr.id());

    gr.initialize(2, 16, 3, bu, gr2, false);
    gr.integrity();
    assertTrue(gr.isInitialized());
    assertEquals(2, gr.sequenceId());
    assertEquals(14, gr.queryEnd());
    assertEquals(12, gr.buildLength());
    assertEquals(expected1, gr.toString());
    assertEquals(bu.bucket(2, 16, 3), gr.bucket());
    assertEquals(bu.bucket(2, 28, 15), gr.bucket());
    assertTrue(gr2 == gr.next());

    gr.write();
    gr.reset();
    gr.integrity();
    assertFalse(gr.isInitialized());
    assertEquals("empty", gr.toString());

    final T  gr3 = getRegion(3, 8, 12);
    gr.initialize(0, 0, 5, bu, gr3, false);
    gr.integrity();
    assertTrue(gr.isInitialized());
    assertEquals(12, gr.buildLength());
    assertEquals(expected2, gr.toString());
    assertEquals(bu.bucket(0, 0, 5), gr.bucket());
    assertEquals(bu.bucket(0, 12, 17), gr.bucket());
    assertTrue(gr3 == gr.next());
  }

  protected void checkNext(final String expected) throws IOException {
    final GapBuckets<T> gb = makeBuckets(8);
    final GapScorer gp = new GappedDistribution(8, 12, GappedDistribution.distrParams(20)).probabilities();

    final T  gr1 = getRegion(1, 8, 12);
    final T  gr2 = getRegion(2, 8, 12);
    final T  gr3 = getRegion(3, 8, 12);
    gr1.initialize(1, 8, 3, gb, gr3, false);
    assertEquals(1, gr1.sequenceId());
    assertEquals(Double.NEGATIVE_INFINITY, gr1.gapScore(gp, gr1));
    assertTrue(gr3 == gr1.next());
    gr2.initialize(1, 24, 19, gb, gr1, false);
    assertTrue(gr1 == gr2.next());
    gr1.setNext(gr1);
    assertTrue(gr1 == gr1.next());
    gr1.setNext(gr2);
    assertTrue(gr2 == gr1.next());
    assertEquals(expected, gr2.toString());
  }

  /** zero delta. */
  protected void checkMerge1(final String expected) throws IOException {
    final GapBuckets<T> gb = makeBuckets(8);
    final GapScorer gp = new GappedDistribution(8, 12, GappedDistribution.distrParams(20)).probabilities();

    final T  gr1 = getRegion(1, 4, 12, gp);
    final T  gr2 = getRegion(2, 4, 12, gp);
    final T  gr3 = getRegion(3, 4, 12, gp);
    final T  gr4 = getRegion(4, 4, 12, gp);
    gr1.initialize(1, 8, 3, gb, gr4, false);
    gr2.initialize(1, 24, 19, gb, gr3, false);
    final long bucket = gr2.bucket();
    assertEquals(bucket, gr1.bucket());
    assertEquals("-0.0359", PROB_FORMAT.format(gr2.gapScore(gp, gr1)));
    gr2.merge(gr1);
    assertTrue(gr1.isInitialized());
    assertTrue(gr2.isInitialized());
    assertEquals(bucket, gr2.bucket());
    assertEquals(28, gr2.buildLength());
    assertEquals(expected, gr2.toString());
  }

  /** +1 delta. */
  protected void checkMerge2a(final String expected) throws IOException {
    final GapBuckets<T> gb = makeBuckets(8);
    final GapScorer gp = new GappedDistribution(8, 12, GappedDistribution.distrParams(20)).probabilities();

    final T  gr1 = getRegion(1, 8, 12);
    final T  gr2 = getRegion(2, 8, 12);
    final T  gr3 = getRegion(3, 8, 12);
    final T  gr4 = getRegion(4, 8, 12);
    gr1.initialize(1, 8, 3, gb, gr3, false);
    gr2.initialize(1, 24, 20, gb, gr4, false);
    final long bucket = gr2.bucket();
    assertEquals(bucket, gr1.bucket());
    assertEquals("-3.8303", PROB_FORMAT.format(gr2.gapScore(gp, gr1)));
    gr2.merge(gr1);
    assertTrue(gr1.isInitialized());
    assertTrue(gr2.isInitialized());
    assertEquals(bucket, gr2.bucket());
    assertEquals(28, gr2.buildLength());
    assertEquals(expected, gr2.toString());
  }

  /** -1 delta. */
  protected void checkMerge3(final String expected) throws IOException {
    final GapBuckets<T> gb = makeBuckets(8);
    final GapScorer gp = new GappedDistribution(8, 12, GappedDistribution.distrParams(20)).probabilities();

    final T  gr1 = getRegion(1, 8, 12);
    final T  gr2 = getRegion(2, 8, 12);
    final T  gr3 = getRegion(3, 8, 12);
    final T  gr4 = getRegion(4, 8, 12);
    gr1.initialize(1, 8, 3, gb, gr3, false);
    gr2.initialize(1, 24, 18, gb, gr4, false);
    final long bucket = gr2.bucket();
    assertEquals(bucket, gr1.bucket());
    assertEquals("-4.0444", PROB_FORMAT.format(gr2.gapScore(gp, gr1)));
    gr2.merge(gr1);
    assertTrue(gr1.isInitialized());
    assertTrue(gr2.isInitialized());
    assertEquals(bucket, gr2.bucket());
    assertEquals(28, gr2.buildLength());
    assertEquals(expected, gr2.toString());
  }

  /** Overlap. */
  protected void checkMerge4(final String expected) throws IOException {
    final GapBuckets<T> gb = makeBuckets(8);
    final GapScorer gp = new GappedDistribution(8, 12, GappedDistribution.distrParams(20)).probabilities();

    final T  gr1 = getRegion(1, 8, 12);
    final T  gr2 = getRegion(2, 8, 12);
    final T  gr3 = getRegion(3, 8, 12);
    final T  gr4 = getRegion(4, 8, 12);
    gr1.initialize(1, 8, 3, gb, gr3, false);
    gr2.initialize(1, 16, 11, gb, gr4, false);
    final long bucket = gr2.bucket();
    assertEquals(bucket, gr1.bucket());
    assertEquals(" 0.0000", PROB_FORMAT.format(gr2.gapScore(gp, gr1)));
    assertEquals(" 0.0000", PROB_FORMAT.format(gr2.gapScore(gr1)));
    gr2.merge(gr1);
    assertTrue(gr1.isInitialized());
    assertTrue(gr2.isInitialized());
    assertEquals(bucket, gr2.bucket());
    assertEquals(20, gr2.buildLength());
    assertEquals(expected, gr2.toString());
  }

  /** Overlap but incorrect. */
  public void testMerge5() throws IOException {
    final GapBuckets<T> gb = makeBuckets(8);
    final GapScorer gp = new GappedDistribution(8, 12, GappedDistribution.distrParams(20)).probabilities();

    final T  gr1 = getRegion(1, 8, 12);
    final T  gr2 = getRegion(2, 8, 12);
    final T  gr3 = getRegion(3, 8, 12);
    final T  gr4 = getRegion(4, 8, 12);
    gr1.initialize(1, 8, 3, gb, gr3, false);
    gr2.initialize(1, 16, 12, gb, gr4, false);
    final long bucket = gr2.bucket();
    assertEquals(bucket, gr1.bucket());
    assertEquals(Double.NEGATIVE_INFINITY, gr2.gapScore(gp, gr1));
  }

  /** Overlap one a subsequence - should merge but not be truncated (bug seen in wild). */
  public void checkMerge6(final String expected) throws IOException {
    final GapBuckets<T> gb = makeBuckets(4);
    final GapScorer gp = new GappedDistribution(4, 12, GappedDistribution.distrParams(20)).probabilities();

    final T  gr1 = getRegion(1, 4, 12);
    final T  gr2 = getRegion(2, 4, 12);
    final T  gr3 = getRegion(3, 4, 12);
    final T  gr4 = getRegion(4, 4, 12);
    gr1.initialize(1, 0, 0, gb, gr3, false);
    gr2.initialize(1, 12, 12, gb, gr4, false);
    final long bucket = gr2.bucket();
    assertEquals(bucket, gr1.bucket());
    assertEquals(0.0, gr2.gapScore(gp, gr1));
    assertEquals(" 0.0000", PROB_FORMAT.format(gr2.gapScore(gr1)));
    gr2.merge(gr1);

    final T  gr5 = getRegion(5, 4, 12);
    final T  gr6 = getRegion(6, 4, 12);
    gr5.initialize(1, 8, 8, gb, gr6, false);
    assertEquals(" 0.0000", PROB_FORMAT.format(gr2.gapScore(gr5)));
    gr2.merge(gr5);
    assertEquals(expected, gr2.toString());
  }

  /** Overlap one a subsequence - should merge but not be truncated (bug seen in wild). */
  public void checkMerge7(final String expected) throws IOException {
    final GapBuckets<T> gb = makeBuckets(4);
    final GapScorer gp = new GappedDistribution(4, 12, GappedDistribution.distrParams(20)).probabilities();

    final T  gr1 = getRegion(1, 4, 12);
    final T  gr2 = getRegion(2, 4, 12);
    final T  gr3 = getRegion(3, 4, 12);
    final T  gr4 = getRegion(4, 4, 12);
    gr1.initialize(1, 0, 0, gb, gr3, false);
    gr2.initialize(1, 12, 12, gb, gr4, false);
    final long bucket = gr2.bucket();
    assertEquals(bucket, gr1.bucket());
    assertEquals(0.0, gr2.gapScore(gp, gr1));
    assertEquals(" 0.0000", PROB_FORMAT.format(gr2.gapScore(gr1)));
    gr2.merge(gr1);

    final T  gr5 = getRegion(5, 4, 12);
    final T  gr6 = getRegion(6, 4, 12);
    gr5.initialize(1, 4, 4, gb, gr6, false);
    assertEquals(" 0.0000", PROB_FORMAT.format(gr2.gapScore(gr5)));
    gr2.merge(gr5);
    assertEquals(expected, gr2.toString());
  }

  /** Special case that gives negative j in probability calculation. */
  public void testGapProbability() throws IOException {
    final GapBuckets<T> gb = makeBuckets(8);
    final GapScorer gp = new GappedDistribution(8, 12, GappedDistribution.distrParams(20)).probabilities();

    final T  gr1 = getRegion(1, 8, 12);
    final T  gr2 = getRegion(2, 8, 12);
    final T  gr3 = getRegion(3, 8, 12);
    final T  gr4 = getRegion(4, 8, 12);
    gr1.initialize(1, 8, 0, gb, gr3, false);
    gr2.initialize(1, 24, 0, gb, gr4, false);
    assertEquals(Double.NEGATIVE_INFINITY, gr2.gapScore(gp, gr1));
  }

  /** Different sequences */
  public void testDifferentSequences() throws IOException {
    final GapBuckets<T> gb = makeBuckets(8);
    final GapScorer gp = new GappedDistribution(8, 12, GappedDistribution.distrParams(20)).probabilities();

    final T  gr1 = getRegion(1, 8, 12);
    final T  gr2 = getRegion(2, 8, 12);
    final T  gr3 = getRegion(3, 8, 12);
    final T  gr4 = getRegion(4, 8, 12);
    gr1.initialize(1, 8, 3, gb, gr3, false);
    gr2.initialize(2, 24, 19, gb, gr4, false);
    final long bucket = gr2.bucket();
    assertNotSame(bucket, gr1.bucket());
    assertEquals(Double.NEGATIVE_INFINITY, gr2.gapScore(gp, gr1));
  }

  /** Big gap */
  public void testBigGap() throws IOException {
    final GapBuckets<T> gb = makeBuckets(8);
    final GapScorer gp = new GappedDistribution(8, 12, GappedDistribution.distrParams(20)).probabilities();

    final T  gr1 = getRegion(1, 8, 12);
    final T  gr2 = getRegion(2, 8, 12);
    final T  gr3 = getRegion(3, 8, 12);
    final T  gr4 = getRegion(4, 8, 12);
    gr1.initialize(1, 8, 3, gb, gr3, false);
    gr2.initialize(1, 48, 43, gb, gr4, false);
    final long bucket = gr2.bucket();
    assertEquals(bucket, gr1.bucket());
    assertEquals(Double.NEGATIVE_INFINITY, gr2.gapScore(gp, gr1));
  }

  public void testBuckets() throws IOException {
    final GapBuckets<T> gb = makeBuckets(8);
    final long low = gb.firstBucket(1);
    final long hi = gb.lastBucket(1);

    final T  gr1 = getRegion(1, 8, 12);
    final T  gr2 = getRegion(2, 8, 12);
    gr1.initialize(1, 8, 3, gb, gr2, false);

    assertEquals(low, gr1.firstBucket());
    assertEquals(hi, gr1.lastBucket());
    final long b = gr1.bucket();
    assertTrue(low <= b && b <= hi);
    for (int i = -5; i < 3; ++i) {
      assertEquals("i=" + i + " b=" + b + " bi=" + gr1.bucket(i), b, gr1.bucket(i));
    }
    assertTrue(b != gr1.bucket(3));
    assertTrue(b != gr1.bucket(-6));
  }

}
