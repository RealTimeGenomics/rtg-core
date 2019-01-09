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

import static com.rtg.util.StringUtils.LS;

import java.io.File;
import java.io.IOException;
import java.io.StringWriter;

import com.rtg.launcher.BuildParams;
import com.rtg.launcher.ISequenceParams;
import com.rtg.launcher.MockReaderParams;
import com.rtg.launcher.MockSequenceParams;
import com.rtg.launcher.ReaderParams;
import com.rtg.mode.BidirectionalFrame;
import com.rtg.mode.SequenceMode;
import com.rtg.reader.MockArraySequencesReader;
import com.rtg.reader.MockSequencesReader;
import com.rtg.reader.SequencesReader;
import com.rtg.util.array.ImmutableIntArray;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.LogRecord;
import com.rtg.util.io.LogStream;

import junit.framework.TestCase;

/**
 */
public class GappedOutputTest extends TestCase {

  private static final String HEADER = "#query-id\tquery-frame\tquery-start\tquery-end\tsubject-id\tsubject-frame\tsubject-start\tsubject-end" + LS;

  protected GappedOutput<GappedRegion> getOutput(final PositionParams params, final Appendable out, final Appendable ambiguousOut, final Appendable unmappedOut) throws IOException {
    final GapBucketsInfo bucketInfo = new GapBucketsInfo(params, 1);
    final PositionWriter writer;
    writer = new RawPositionWriter(out, unmappedOut);
    final GappedRegionFactory<GappedRegion> regionFactory = new GappedRegionFactory<GappedRegion>() {
      @Override
      public GappedRegion region(int id, GapScorer distribution, PositionParams params, ImmutableIntArray buildLengths) {
        return new GappedRegion(id, params.build(), distribution);
      }
    };
    return new GappedOutput<>(params, regionFactory, new GappedDistribution(params).probabilities(), null, out, writer, bucketInfo);
  }

  protected static PositionParams makeParams(final int[] seqLengths, final int stepSize, final int windowSize, final int gapSize) throws IOException {
    return makeParams(seqLengths, stepSize, windowSize, gapSize, 3/*threshold*/);
  }

  protected static PositionParams makeParams(final int[] seqLengths, final int stepSize, final int windowSize, final int gapSize, final int threshold) throws IOException {
    final SequencesReader reader = new MockArraySequencesReader(SequenceMode.UNIDIRECTIONAL.type(), seqLengths);
    final ReaderParams srp = new MockReaderParams(reader);
    final ISequenceParams subjectParams = new MockSequenceParams(srp, SequenceMode.UNIDIRECTIONAL, 0, reader.numberSequences());
    final BuildParams buildParams = BuildParams.builder().windowSize(windowSize).stepSize(stepSize).sequences(subjectParams).create();

    final SequencesReader reader2 = new MockSequencesReader(SequenceMode.BIDIRECTIONAL.codeType(), 50);
    final ReaderParams qrp = new MockReaderParams(reader2);
    final ISequenceParams queryParams = new MockSequenceParams(qrp, SequenceMode.BIDIRECTIONAL, 0, 49);
    final BuildParams searchParams = BuildParams.builder().windowSize(windowSize).stepSize(1).sequences(queryParams).create();
    final PositionDistributionParams distr = new PositionDistributionParams(0.001, 0.009, gapSize, 0);
    final PositionOutputParams outParams = new PositionOutputParams(new File(""), OutputFormatType.SEGMENT, distr, null, false, 0);
    return PositionParams.builder().hashCountThreshold(threshold).buildParams(buildParams).searchParams(searchParams).outputParams(outParams).create();
  }

  @Override
  public void setUp() {
    Diagnostic.setLogStream();
  }

  /**
   * Simple succesive words.
   * @throws IOException
   */
  public void test1() throws IOException {
    final PositionParams params = makeParams(new int[] {12}, 4/*step*/, 4/*word*/, 8/*gap*/);
    final Appendable out = new StringWriter();
    final Appendable unmappedOut = new StringWriter();
    final AbstractPositionOutput so = getOutput(params, out, unmappedOut, unmappedOut);
    so.globalIntegrity();
    final String emptyExpected = ""
      + "GappedOutput wordSize=4 maxGap=8 minDelta=-1 maxDelta=1 score=" + 0.0 + LS
      + "[0]0" + LS
      + "[1]0" + LS
      + "[2]0" + LS
      + "[3]0" + LS
      + "[4]0" + LS
      + "[5]0" + LS
      + "[6]0" + LS
      + "[7]0" + LS
      + "[8]0" + LS
      ;
    assertEquals(emptyExpected, so.toString());
    so.nextSequence(42, 0);
    so.nextQuery(BidirectionalFrame.FORWARD, 42);

    posnHit(so, 0, 0);
    posnHit(so, 4, 4);
    posnHit(so, 8, 8);
    posnHit(so, 12, 12);
    posnHit(so, 16, 16);

    final String fullExpected = ""
      + "GappedOutput wordSize=4 maxGap=8 minDelta=-1 maxDelta=1 score=" + 0.0 + LS
      + "[0]0" + LS
      + "[1]0" + LS
      + "[2]0" + LS
      + "[3]1 4[12..15]0:[12..15]^0->5" + LS
      + "[4]0" + LS
      + "[5]0" + LS
      + "[6]0" + LS
      + "[7]1 5[16..19]0:[16..19]^0->3" + LS
      + "[8]1 3[0..11]0:[0..11]^0->4" + LS
      ;
    assertEquals(fullExpected, so.toString());

    so.endQuery();
    so.globalIntegrity();
    so.endQuerySequence();
    so.globalIntegrity();
    so.endAll();
    so.globalIntegrity();
    final String expected =  HEADER + "42\tF\t1\t20\t0\t\t1\t20" + LS;
    assertEquals(expected, out.toString());
    assertEquals(20.0, so.score());
  }

  /**
   * Single gap.
   * @throws IOException
   */
  public void test2() throws IOException {
    final PositionParams params = makeParams(new int[] {12}, 4/*step*/, 4/*word*/, 8/*gap*/);
    final Appendable out = new StringWriter();
    final Appendable unmappedOut = new StringWriter();
    final AbstractPositionOutput so = getOutput(params, out, unmappedOut, unmappedOut);
    so.globalIntegrity();
    so.nextSequence(42, 0);
    so.nextQuery(BidirectionalFrame.FORWARD, 42);

    posnHit(so, 0, 0);
    posnHit(so, 4, 4);
    posnHit(so, 12, 12);
    posnHit(so, 16, 16);
    //System.err.println(so.toString());
    so.endQuery();
    so.globalIntegrity();
    so.endQuerySequence();
    so.globalIntegrity();
    so.endAll();
    so.globalIntegrity();
    final String expected =  HEADER + "42\tF\t1\t20\t0\t\t1\t20" + LS;
    assertEquals(expected, out.toString());
    assertEquals(20.0, so.score());
  }

  /**
   * Two gaps.
   * @throws IOException
   */
  public void test3() throws IOException {
    final PositionParams params = makeParams(new int[] {12}, 4/*step*/, 4/*word*/, 8/*gap*/);
    final Appendable out = new StringWriter();
    final Appendable unmappedOut = new StringWriter();
    final AbstractPositionOutput so = getOutput(params, out, unmappedOut, unmappedOut);
    so.globalIntegrity();
    so.nextSequence(42, 0);
    so.nextQuery(BidirectionalFrame.FORWARD, 42);

    posnHit(so, 0, 0);
    posnHit(so, 8, 8);
    posnHit(so, 16, 16);

    so.endQuery();
    so.globalIntegrity();
    so.endQuerySequence();
    so.globalIntegrity();
    so.endAll();
    so.globalIntegrity();
    final String expected =  HEADER + "42\tF\t1\t20\t0\t\t1\t20" + LS;
    assertEquals(expected, out.toString());
    assertEquals(20.0, so.score());
  }

  /**
   * Single insertion.
   * @throws IOException
   */
  public void test4() throws IOException {
    final PositionParams params = makeParams(new int[] {12}, 4/*step*/, 4/*word*/, 8/*gap*/);
    final Appendable out = new StringWriter();
    final Appendable unmappedOut = new StringWriter();
    final AbstractPositionOutput so = getOutput(params, out, unmappedOut, unmappedOut);
    so.globalIntegrity();
    so.nextSequence(42, 0);
    so.nextQuery(BidirectionalFrame.FORWARD, 42);

    posnHit(so, 0, 0);
    posnHit(so, 4, 4);
    posnHit(so, 12, 11);
    posnHit(so, 16, 15);

    so.endQuery();
    so.globalIntegrity();
    so.endQuerySequence();
    so.globalIntegrity();
    so.endAll();
    so.globalIntegrity();
    final String expected =  HEADER + "42\tF\t1\t19\t0\t\t1\t20" + LS;
    assertEquals(expected, out.toString());
    assertEquals(20.0, so.score());
  }


  /**
   * Deletion and substitution - longer gap will be missed.
   * @throws IOException
   */
  public void test5() throws IOException {
    final PositionParams params = makeParams(new int[] {12}, 4/*step*/, 4/*word*/, 8/*gap*/);
    final Appendable out = new StringWriter();
    final Appendable unmappedOut = new StringWriter();
    final AbstractPositionOutput so = getOutput(params, out, unmappedOut, unmappedOut);
    so.globalIntegrity();
    so.nextSequence(42, 0);
    so.nextQuery(BidirectionalFrame.FORWARD, 42);

    posnHit(so, 0, 0);
    posnHit(so, 4, 4);
    posnHit(so, 12, 13);
    posnHit(so, 16, 17);

    so.endQuery();
    so.globalIntegrity();
    so.endQuerySequence();
    so.globalIntegrity();
    so.endAll();
    so.globalIntegrity();
    final String expected =  HEADER + "42\tF\t1\t8\t0\t\t1\t8" + LS + "42\tF\t14\t21\t0\t\t13\t20" + LS;
    assertEquals(expected, out.toString());
    assertEquals(16.0, so.score());
  }

  /**
   * Deletion and substitution - longer gap.
   * @throws IOException
   */
  public void test6() throws IOException {
    final PositionParams params = makeParams(new int[] {12}, 4/*step*/, 4/*word*/, 12/*gap*/);
    final Appendable out = new StringWriter();
    final Appendable unmappedOut = new StringWriter();
    final AbstractPositionOutput so = getOutput(params, out, unmappedOut, unmappedOut);
    so.globalIntegrity();
    so.nextSequence(42, 0);
    so.nextQuery(BidirectionalFrame.FORWARD, 42);

    posnHit(so, 0, 0);
    posnHit(so, 4, 4);
    posnHit(so, 12, 13);
    posnHit(so, 16, 17);

    so.endQuery();
    so.endQuerySequence();
    so.endAll();
    so.globalIntegrity();
    final String expected =  HEADER + "42\tF\t1\t21\t0\t\t1\t20" + LS;
    assertEquals(20.0, so.score());
    assertEquals(expected, out.toString());
  }

  /**
   * Overlap.
   * @throws IOException
   */
  public void test7() throws IOException {
    final PositionParams params = makeParams(new int[] {12}, 2/*step*/, 4/*word*/, 12/*gap*/);
    final Appendable out = new StringWriter();
    final Appendable unmappedOut = new StringWriter();
    final AbstractPositionOutput so = getOutput(params, out, unmappedOut, unmappedOut);
    so.globalIntegrity();
    so.nextSequence(42, 0);
    so.nextQuery(BidirectionalFrame.FORWARD, 42);

    posnHit(so, 0, 0);
    posnHit(so, 2, 2);

    so.endQuery();
    so.endQuerySequence();
    so.endAll();
    so.globalIntegrity();
    final String expected =  HEADER + "42\tF\t1\t6\t0\t\t1\t6" + LS;
    assertEquals(expected, out.toString());
    assertEquals(6.0, so.score());
  }

  /**
   * Single word hit at position 1 - suggested by jumble.
   * @throws IOException
   */
  public void test8() throws IOException {
    final PositionParams params = makeParams(new int[] {12}, 2/*step*/, 4/*word*/, 12/*gap*/);
    final Appendable out = new StringWriter();
    final Appendable unmappedOut = new StringWriter();
    final AbstractPositionOutput so = getOutput(params, out, unmappedOut, unmappedOut);
    so.globalIntegrity();
    so.nextSequence(42, 0);
    so.nextQuery(BidirectionalFrame.FORWARD, 42);

    posnHit(so, 0, 1);

    so.endQuery();
    so.endQuerySequence();
    so.endAll();
    so.globalIntegrity();
    final String expected =  HEADER + "42\tF\t2\t5\t0\t\t1\t4" + LS;
    assertEquals(expected, out.toString());
    assertEquals(4.0, so.score());
  }

  /**
   * Single word hit at position 1, followed by another hit - suggested by jumble.
   * @throws IOException
   */
  public void test9() throws IOException {
    final PositionParams params = makeParams(new int[] {12}, 2/*step*/, 4/*word*/, 12/*gap*/);
    final Appendable out = new StringWriter();
    final Appendable unmappedOut = new StringWriter();
    final AbstractPositionOutput so = getOutput(params, out, unmappedOut, unmappedOut);
    so.globalIntegrity();
    so.nextSequence(42, 0);
    so.nextQuery(BidirectionalFrame.FORWARD, 42);

    posnHit(so, 0, 1);
    posnHit(so, 4, 5);

    so.endQuery();
    so.endQuerySequence();
    so.endAll();
    so.globalIntegrity();
    final String expected =  HEADER + "42\tF\t2\t9\t0\t\t1\t8" + LS;
    assertEquals(expected, out.toString());
    assertEquals(8.0, so.score());
  }

  /**
   * From case that caused problems in Position unit tests.
   * Ensures that multiple bucket chains are compared which gives some interesting cases when computing
   * the gap probabilities.
   * @throws IOException
   */
  public void test10() throws IOException {
    final PositionParams params = makeParams(new int[] {12}, 4/*step*/, 4/*word*/, 12/*gap*/);
    final Appendable out = new StringWriter();
    final Appendable unmappedOut = new StringWriter();
    final AbstractPositionOutput so = getOutput(params, out, unmappedOut, unmappedOut);
    so.globalIntegrity();
    so.nextSequence(42, 0);
    so.nextQuery(BidirectionalFrame.FORWARD, 0);

    so.setPosition(0);
    so.hit(0, 0);
    so.hit(0, 4);
    so.hit(0, 8);
    so.globalIntegrity();
    so.endPosition();
    so.globalIntegrity();

    so.setPosition(9);
    so.hit(0, 0);
    so.hit(0, 4);
    so.hit(0, 8);
    so.globalIntegrity();
    so.endPosition();
    so.globalIntegrity();

    so.endQuery();
    so.endQuerySequence();
    so.endAll();
    so.globalIntegrity();
    final String expected = HEADER
      + "0\tF\t1\t4\t0\t\t5\t8" + LS
      + "0\tF\t1\t4\t0\t\t9\t12" + LS
      + "0\tF\t10\t13\t0\t\t1\t4" + LS
      + "0\tF\t10\t13\t0\t\t5\t8" + LS
      + "0\tF\t1\t13\t0\t\t1\t12" + LS
      ;
    assertEquals(expected, out.toString());
    assertEquals(28.0, so.score());
  }

  public void testLog() throws IOException {
    final LogStream logStream = new LogRecord();
    Diagnostic.setLogStream(logStream);
    final PositionParams params = makeParams(new int[] {12}, 2/*step*/, 4/*word*/, 12/*gap*/);
    final Appendable out = new StringWriter();
    final Appendable unmappedOut = new StringWriter();
    final AbstractPositionOutput so = getOutput(params, out, unmappedOut, unmappedOut);
    so.globalIntegrity();
    so.nextSequence(42, 0);
    so.nextQuery(BidirectionalFrame.FORWARD, 42);

    posnHit(so, 0, 1);

    so.endQuery();
    so.endQuerySequence();
    so.endAll();
    so.globalIntegrity();
    final String expected =  HEADER + "42\tF\t2\t5\t0\t\t1\t4" + LS;
    assertEquals(expected, out.toString());
    final String log = logStream.toString();
    //System.err.println(log);
    assertTrue(log.contains("GappedOutput minDelta=-4 maxDelta=4"));
  }

  private void posnHit(final AbstractPositionOutput so, final int bp, final int qp) throws IOException {
    so.setPosition(qp);
    so.hit(0, bp);
    so.globalIntegrity();
    so.endPosition();
    so.globalIntegrity();
  }

  /**
   * Big complicated case - will fill the regions array.
   * All simple adjacent words.
   * @throws IOException
   */
  public void testBig() throws IOException {
    final PositionParams params = makeParams(new int[] {67, 145, 220}, 4/*step*/, 4/*word*/, 8/*gap*/);
    checkBig(params);
  }


  //exercise case when the bucket array is smaller than the range of bucket values (so there will be collisions in buckets)
  public void testBigBucketHash() throws IOException {
    final PositionParams params = makeParams(new int[] {67, 145, 220}, 1/*step*/, 4/*word*/, 8/*gap*/, 3/*threshold*/);
    checkBig(params);
  }

  private void checkBig(final PositionParams params) throws IOException {
    final Appendable out = new StringWriter();
    final Appendable unmappedOut = new StringWriter();
    final AbstractPositionOutput so = getOutput(params, out, unmappedOut, unmappedOut);
    so.globalIntegrity();
    so.nextSequence(42, 0);
    so.nextQuery(BidirectionalFrame.REVERSE, 43);

    for (int i = 0; i < 12; i += 4) {
      so.setPosition(i);
      so.hit(0, i + 4);
      so.hit(1, i);
      so.hit(2, i + 8);
      so.globalIntegrity();
      so.endPosition();
      so.globalIntegrity();

      so.setPosition(i + 1);
      so.hit(0, i + 20);
      so.hit(1, i + 12);
      so.hit(2, i + 24);
      so.globalIntegrity();
      so.endPosition();
      so.globalIntegrity();

      so.setPosition(i + 2);
      so.hit(0, i + 36);
      so.hit(1, i + 48);
      so.hit(2, i + 72);
      so.globalIntegrity();
      so.endPosition();
      so.globalIntegrity();

      so.setPosition(i + 3);
      so.hit(0, i + 16);
      so.hit(1, i + 56);
      so.hit(2, i + 44);
      so.globalIntegrity();
      so.endPosition();
      so.globalIntegrity();
    }

    so.endQuery();
    so.endQuerySequence();
    so.endAll();
    so.globalIntegrity();
    final String expected =  HEADER
      + "43\tR\t1\t12\t0\t\t5\t16" + LS
      + "43\tR\t1\t12\t1\t\t1\t12" + LS
      + "43\tR\t1\t12\t2\t\t9\t20" + LS
      + "43\tR\t2\t13\t0\t\t21\t32" + LS
      + "43\tR\t2\t13\t1\t\t13\t24" + LS
      + "43\tR\t2\t13\t2\t\t25\t36" + LS
      + "43\tR\t3\t14\t0\t\t37\t48" + LS
      + "43\tR\t3\t14\t1\t\t49\t60" + LS
      + "43\tR\t3\t14\t2\t\t73\t84" + LS
      + "43\tR\t4\t15\t0\t\t17\t28" + LS
      + "43\tR\t4\t15\t1\t\t57\t68" + LS
      + "43\tR\t4\t15\t2\t\t45\t56" + LS
      ;
    assertEquals(expected, out.toString());
    assertEquals(144.0, so.score());
  }
}
