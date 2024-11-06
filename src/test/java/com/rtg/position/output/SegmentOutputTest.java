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
package com.rtg.position.output;

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
import com.rtg.util.StringUtils;

import junit.framework.TestCase;

/**
 */
public class SegmentOutputTest extends TestCase {

  /**
   * Header
   */
  public static final String HEADER = "#query-id\t"
    + "query-frame\t"
    + "query-start\t"
    + "subject-id\t"
    + "subject-frame\t"
    + "subject-start\t"
    + "length" + com.rtg.util.StringUtils.LS;

  protected AbstractPositionOutput getOutput(final PositionParams params, final Appendable out) {
    return new SegmentOutput(params, out);
  }

  protected static PositionParams makeParams(final int threshold, final int windowSize, final int stepSize) throws IOException {
    final ReaderParams srp = new MockReaderParams(0, 0, SequenceMode.UNIDIRECTIONAL.codeType());
    final ISequenceParams subjectParams = new MockSequenceParams(srp, SequenceMode.UNIDIRECTIONAL, 0, 0);
    final BuildParams buildParams = BuildParams.builder().windowSize(windowSize).stepSize(stepSize).sequences(subjectParams).create();

    final ReaderParams qrp = new MockReaderParams(0, 0, SequenceMode.BIDIRECTIONAL.codeType());
    final ISequenceParams queryParams = new MockSequenceParams(qrp, SequenceMode.BIDIRECTIONAL, 0, 0);
    final BuildParams searchParams = BuildParams.builder().windowSize(windowSize).stepSize(1).sequences(queryParams).create();
    final PositionOutputParams outParams = new PositionOutputParams(new File(""), OutputFormatType.SEGMENT, null, null, false, 10);
    return PositionParams.builder().hashCountThreshold(threshold).buildParams(buildParams).searchParams(searchParams).outputParams(outParams).create();
  }

  /**
   * Special case when windowSize == 2.
   * @throws IOException
   */
  public void test1a() throws IOException {
    final PositionParams params = makeParams(3, 4, 4);
    final Appendable out = new StringWriter();
    final AbstractPositionOutput so = getOutput(params, out);
    so.globalIntegrity();
    assertEquals("SegmentOutput", so.toString());
    so.nextSequence(42, 0);
    so.nextQuery(BidirectionalFrame.FORWARD, 42);

    so.setPosition(3);
    so.hit(7, 0);
    so.hit(7, 1);
    so.hit(7, 2);
    so.globalIntegrity();
    so.endPosition();
    so.globalIntegrity();

    so.setPosition(5);
    so.hit(9, 0);
    so.hit(9, 1);
    so.hit(9, 2);
    so.globalIntegrity();
    so.endPosition();
    so.globalIntegrity();

    so.setPosition(6);
    so.hit(10, 0);
    so.hit(10, 1);
    so.hit(10, 2);
    so.globalIntegrity();
    so.endPosition();
    so.globalIntegrity();

    so.setPosition(11);
    so.hit(5, 0);
    so.hit(5, 1);
    so.hit(5, 2);
    so.globalIntegrity();
    so.endPosition();
    so.globalIntegrity();

    so.setPosition(23);
    so.hit(5, 0);
    so.hit(5, 1);
    so.hit(5, 2);
    so.globalIntegrity();
    so.endPosition();
    so.globalIntegrity();

    so.endQuery();
    so.globalIntegrity();
    so.endQuerySequence();
    so.globalIntegrity();
    so.endAll();
    so.globalIntegrity();
    final String expected = ""
      + "42\tF\t4\t7\t\t1\t4" + StringUtils.LS
      + "42\tF\t4\t7\t\t2\t4" + StringUtils.LS
      + "42\tF\t4\t7\t\t3\t4" + StringUtils.LS
      + "42\tF\t6\t9\t\t1\t4" + StringUtils.LS
      + "42\tF\t6\t9\t\t2\t4" + StringUtils.LS
      + "42\tF\t6\t9\t\t3\t4" + StringUtils.LS
      + "42\tF\t7\t10\t\t1\t4" + StringUtils.LS
      + "42\tF\t7\t10\t\t2\t4" + StringUtils.LS
      + "42\tF\t7\t10\t\t3\t4" + StringUtils.LS
      + "42\tF\t12\t5\t\t1\t4" + StringUtils.LS
      + "42\tF\t12\t5\t\t2\t4" + StringUtils.LS
      + "42\tF\t12\t5\t\t3\t4" + StringUtils.LS
      + "42\tF\t24\t5\t\t1\t4" + StringUtils.LS
      + "42\tF\t24\t5\t\t2\t4" + StringUtils.LS
      + "42\tF\t24\t5\t\t3\t4" + StringUtils.LS
      ;
    assertEquals(HEADER + expected, out.toString());
    assertEquals(60.0, so.score());
  }

  /**
   * Small gap and big gap in setPosition.
   * @throws IOException
   */
  public void test1b() throws IOException {
    final PositionParams params = makeParams(2, 3, 3);
    final Appendable out = new StringWriter();
    final AbstractPositionOutput so = getOutput(params, out);
    so.globalIntegrity();
    assertEquals("SegmentOutput", so.toString());
    so.nextSequence(42, 0);
    so.nextQuery(BidirectionalFrame.FORWARD, 42);

    //fill adjacent collections
    so.setPosition(3);
    so.hit(7, 1);
    so.globalIntegrity();
    so.endPosition();
    so.globalIntegrity();

    so.setPosition(4);
    so.hit(1, 11);
    so.hit(2, 11);
    so.globalIntegrity();
    so.endPosition();
    so.globalIntegrity();

    so.setPosition(5);
    so.hit(9, 1);
    so.globalIntegrity();
    so.endPosition();
    so.globalIntegrity();

    //one step jump - see 4
    so.setPosition(7);
    so.hit(1, 14);
    so.globalIntegrity();
    so.endPosition();
    so.globalIntegrity();

    //4 (one more than word size) jump
    so.setPosition(11);
    so.hit(5, 1);
    so.globalIntegrity();
    so.endPosition();
    so.globalIntegrity();

    so.endQuery();
    so.globalIntegrity();
    so.endQuerySequence();
    so.globalIntegrity();
    so.endAll();
    so.globalIntegrity();
    final String expected = ""
      + "42\tF\t4\t7\t\t2\t3" + StringUtils.LS
      + "42\tF\t5\t2\t\t12\t3" + StringUtils.LS
      + "42\tF\t6\t9\t\t2\t3" + StringUtils.LS
      + "42\tF\t5\t1\t\t12\t6" + StringUtils.LS
      + "42\tF\t12\t5\t\t2\t3" + StringUtils.LS
      ;
    assertEquals(HEADER + expected, out.toString());
    assertEquals(18.0, so.score());
  }


  /**
   * Same as test1b except that stepsize = 2
   * Small gap and big gap in setPosition.
   * @throws IOException
   */
  public void test1c() throws IOException {
    final PositionParams params = makeParams(2, 3, 2);
    final Appendable out = new StringWriter();
    final AbstractPositionOutput so = getOutput(params, out);
    so.globalIntegrity();
    assertEquals("SegmentOutput", so.toString());
    so.nextSequence(42, 0);
    so.nextQuery(BidirectionalFrame.FORWARD, 42);

    //fill adjacent collections
    so.setPosition(3);
    so.hit(7, 1);
    so.globalIntegrity();
    so.endPosition();
    so.globalIntegrity();

    so.setPosition(4);
    so.hit(1, 11);
    so.hit(2, 11);
    so.globalIntegrity();
    so.endPosition();
    so.globalIntegrity();

    so.setPosition(5);
    so.hit(9, 1);
    so.globalIntegrity();
    so.endPosition();
    so.globalIntegrity();

    //one step jump - see 4
    so.setPosition(6);
    so.hit(1, 13);
    so.globalIntegrity();
    so.endPosition();
    so.globalIntegrity();

    //3 (one more than word size) jump
    so.setPosition(9);
    so.hit(5, 1);
    so.globalIntegrity();
    so.endPosition();
    so.globalIntegrity();

    so.endQuery();
    so.globalIntegrity();
    so.endQuerySequence();
    so.globalIntegrity();
    so.endAll();
    so.globalIntegrity();
    final String expected = ""
      + "42\tF\t4\t7\t\t2\t3" + StringUtils.LS
      + "42\tF\t5\t2\t\t12\t3" + StringUtils.LS
      + "42\tF\t6\t9\t\t2\t3" + StringUtils.LS
      + "42\tF\t5\t1\t\t12\t5" + StringUtils.LS
      + "42\tF\t10\t5\t\t2\t3" + StringUtils.LS
      ;
    assertEquals(HEADER + expected, out.toString());
    assertEquals(17.0, so.score());
  }

  /**
   * Same as test1b except that stepsize = 2
   * Small gap and big gap in setPosition.
   * @throws IOException
   */
  public void test1d() throws IOException {
    final PositionParams params = makeParams(2, 3, 2);
    final Appendable out = new StringWriter();
    final AbstractPositionOutput so = getOutput(params, out);
    so.globalIntegrity();
    assertEquals("SegmentOutput", so.toString());
    so.nextSequence(42, 0);
    so.nextQuery(BidirectionalFrame.FORWARD, 42);

    //fill adjacent collections
    so.setPosition(3);
    so.hit(7, 1);
    so.hit(8, 1);
    so.globalIntegrity();
    so.endPosition();
    so.globalIntegrity();

    so.setPosition(4);
    so.hit(1, 11);
    so.hit(2, 11);
    so.globalIntegrity();
    so.endPosition();
    so.globalIntegrity();

    //3 (one more than word size) jump
    so.setPosition(7);
    so.hit(7, 3);
    so.globalIntegrity();
    so.endPosition();
    so.globalIntegrity();

    //the critical one that will blow up if setposition isnt just right
    //(see 4)
    so.setPosition(8);
    so.hit(1, 13);
    so.globalIntegrity();
    so.endPosition();
    so.globalIntegrity();

    so.endQuery();
    so.globalIntegrity();
    so.endQuerySequence();
    so.globalIntegrity();
    so.endAll();
    so.globalIntegrity();
    final String expected = ""
      + "42\tF\t4\t7\t\t2\t3" + StringUtils.LS
      + "42\tF\t4\t8\t\t2\t3" + StringUtils.LS
      + "42\tF\t5\t1\t\t12\t3" + StringUtils.LS
      + "42\tF\t5\t2\t\t12\t3" + StringUtils.LS
      + "42\tF\t8\t7\t\t4\t3" + StringUtils.LS
      + "42\tF\t9\t1\t\t14\t3" + StringUtils.LS
      ;
    assertEquals(HEADER + expected, out.toString());
    assertEquals(18.0, so.score());
  }


  /**
   * Three succesive merges (inspired by a bug in calculating the start of the build position on output).
   * @throws IOException
   */
  public void test1e() throws IOException {
    final PositionParams params = makeParams(2, 3, 2);
    final Appendable out = new StringWriter();
    final AbstractPositionOutput so = getOutput(params, out);
    so.globalIntegrity();
    assertEquals("SegmentOutput", so.toString());
    so.nextSequence(42, 0);
    so.nextQuery(BidirectionalFrame.FORWARD, 42);

    //fill adjacent collections
    so.setPosition(3);
    so.hit(7, 1);
    so.globalIntegrity();
    so.endPosition();
    so.globalIntegrity();

    so.setPosition(5);
    so.hit(7, 3);
    so.globalIntegrity();
    so.endPosition();
    so.globalIntegrity();

    //3 (one more than word size) jump
    so.setPosition(7);
    so.hit(7, 5);
    so.globalIntegrity();
    so.endPosition();
    so.globalIntegrity();

    so.endQuery();
    so.globalIntegrity();
    so.endQuerySequence();
    so.globalIntegrity();
    so.endAll();
    so.globalIntegrity();
    final String expected = ""
      + "42\tF\t4\t7\t\t2\t7" + StringUtils.LS
      ;
    assertEquals(HEADER + expected, out.toString());
    assertEquals(7.0, so.score());
  }

  /**
   * Test case when the free list is exhausted briefly.
   * Fill each collection and at endPosition generate new hits before freeing all the old ones.
   * @throws IOException
   */
  public void testFreeList() throws IOException {
    final PositionParams params = makeParams(3, 4, 4);
    final Appendable out = new StringWriter();
    final Appendable unmappedOut = new StringWriter();
    final AbstractPositionOutput so = getOutput(params, out);
    so.globalIntegrity();
    so.nextSequence(42, 0);
    so.nextQuery(BidirectionalFrame.FORWARD, 42);

    so.setPosition(3);
    so.hit(7, 0);
    so.hit(7, 1);
    so.hit(7, 2);
    so.globalIntegrity();
    so.endPosition();
    so.globalIntegrity();

    so.setPosition(4);
    so.hit(8, 0);
    so.hit(8, 1);
    so.hit(8, 2);
    so.globalIntegrity();
    so.endPosition();
    so.globalIntegrity();

    so.setPosition(5);
    so.hit(9, 0);
    so.hit(9, 1);
    so.hit(9, 2);
    so.globalIntegrity();
    so.endPosition();
    so.globalIntegrity();

    so.setPosition(6);
    so.hit(10, 0);
    so.hit(10, 1);
    so.hit(10, 2);
    so.globalIntegrity();
    so.endPosition();
    so.globalIntegrity();

    //should be full
    //force new values less than the ones added for position 3
    so.setPosition(7);
    so.hit(5, 0); //5 < 7
    so.hit(5, 1);
    so.hit(5, 2);
    so.globalIntegrity();
    so.endPosition();
    so.globalIntegrity();


    so.endQuery();
    so.globalIntegrity();
    so.endQuerySequence();
    so.globalIntegrity();
    so.endAll();
    so.globalIntegrity();
    final String expected = ""
      + "42\tF\t4\t7\t\t1\t4" + StringUtils.LS
      + "42\tF\t4\t7\t\t2\t4" + StringUtils.LS
      + "42\tF\t4\t7\t\t3\t4" + StringUtils.LS
      + "42\tF\t5\t8\t\t1\t4" + StringUtils.LS
      + "42\tF\t5\t8\t\t2\t4" + StringUtils.LS
      + "42\tF\t5\t8\t\t3\t4" + StringUtils.LS
      + "42\tF\t6\t9\t\t1\t4" + StringUtils.LS
      + "42\tF\t6\t9\t\t2\t4" + StringUtils.LS
      + "42\tF\t6\t9\t\t3\t4" + StringUtils.LS
      + "42\tF\t7\t10\t\t1\t4" + StringUtils.LS
      + "42\tF\t7\t10\t\t2\t4" + StringUtils.LS
      + "42\tF\t7\t10\t\t3\t4" + StringUtils.LS
      + "42\tF\t8\t5\t\t1\t4" + StringUtils.LS
      + "42\tF\t8\t5\t\t2\t4" + StringUtils.LS
      + "42\tF\t8\t5\t\t3\t4" + StringUtils.LS
      ;
    assertEquals(HEADER + expected, out.toString());
    assertEquals(60.0, so.score());
    assertEquals("", unmappedOut.toString());
  }

  /**
   * Test case when no hits but other things happen.
   * @throws IOException
   */
  public void testNoHits() throws IOException {
    final PositionParams params = makeParams(3, 4, 4);
    final Appendable out = new StringWriter();
    final AbstractPositionOutput so = getOutput(params, out);
    so.globalIntegrity();
    so.nextSequence(42, 0);
    so.nextQuery(BidirectionalFrame.FORWARD, 42);

    so.setPosition(3);
    so.globalIntegrity();
    so.endPosition();
    so.globalIntegrity();

    so.setPosition(4);
    so.globalIntegrity();
    so.endPosition();
    so.globalIntegrity();

    so.setPosition(5);
    so.globalIntegrity();
    so.endPosition();
    so.globalIntegrity();

    so.setPosition(6);
    so.globalIntegrity();
    so.endPosition();
    so.globalIntegrity();

    so.setPosition(7);
    so.globalIntegrity();
    so.endPosition();
    so.globalIntegrity();


    so.endQuery();
    so.globalIntegrity();
    so.endQuerySequence();
    so.globalIntegrity();
    so.endAll();
    so.globalIntegrity();
    assertEquals("", out.toString());
    assertEquals(0.0, so.score());
  }

}
