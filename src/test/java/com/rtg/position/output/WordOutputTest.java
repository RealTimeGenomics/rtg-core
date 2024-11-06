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

import java.io.IOException;
import java.io.StringWriter;

import com.rtg.mode.BidirectionalFrame;
import com.rtg.mode.UnidirectionalFrame;
import com.rtg.util.StringUtils;

import junit.framework.TestCase;



/**
 */
public class WordOutputTest extends TestCase {

  private static final String EXPECTED1 = ""
    + "0\tF\t1\t42\t\t14" + StringUtils.LS
    + "0\tF\t1\t17\t\t20" + StringUtils.LS
    + "0\tF\t54\t17\t\t24" + StringUtils.LS
    + "0\tR\t1\t43\t\t13" + StringUtils.LS
    ;

  /**
   * HEader
   */
  public static final String HEADER = "#query-id\tquery-frame\tquery-start\tsubject-id\tsubject-frame\tsubject-start" + StringUtils.LS;
  public final void test1() throws IOException {
    final Appendable out = new StringWriter();
    final Appendable unmappedOut = new StringWriter();
    final PositionOutput wo
    =  new WordOutput(UnidirectionalFrame.values(), out, unmappedOut);
    wo.nextSequence(0, 0);
    wo.nextQuery(BidirectionalFrame.FORWARD, 0);
    wo.setPosition(0);
    //the actual values of the ints here are not important (well they need to be different)
    wo.hit(42, 13);
    wo.hit(17, 19);
    wo.endPosition();

    wo.setPosition(51);
    wo.endPosition();

    wo.setPosition(53);
    wo.hit(17, 23);
    wo.endPosition();
    wo.endQuery();

    wo.nextQuery(BidirectionalFrame.REVERSE, 0);
    wo.setPosition(0);
    wo.hit(43, 12);
    wo.endPosition();
    wo.endQuery();
    wo.endQuerySequence();
    wo.endAll();
    assertEquals(HEADER + EXPECTED1, out.toString());
    assertEquals(4.0, wo.score());
    assertEquals("", unmappedOut.toString());
  }


  private static final String EXPECTED2 = ""
    + "0\tF\t1\t21\tF\t14" + StringUtils.LS
    + "0\tF\t1\t8\tR\t20" + StringUtils.LS
    + "0\tF\t54\t8\tR\t24" + StringUtils.LS
    + "0\tR\t1\t21\tR\t13" + StringUtils.LS
    ;

  /**
   * Bidirectional (tests the case when number frames > 1).
   * start position > 0.
   */
  public final void test2() throws IOException {
    final Appendable out = new StringWriter();
    final Appendable unmappedOut = new StringWriter();
    final PositionOutput wo =  new WordOutput(BidirectionalFrame.values(), out, unmappedOut);
    wo.nextSequence(0, 0);
    wo.nextQuery(BidirectionalFrame.FORWARD, 0);
    wo.setPosition(0);
    //the actual values of the ints here are not important (well they need to be different)
    wo.hit(42, 13);
    wo.hit(17, 19);
    wo.endPosition();

    wo.setPosition(51);
    wo.endPosition();

    wo.setPosition(53);
    wo.hit(17, 23);
    wo.endPosition();
    wo.endQuery();

    wo.nextQuery(BidirectionalFrame.REVERSE, 0);
    wo.setPosition(0);
    wo.hit(43, 12);
    wo.endPosition();
    wo.endQuery();
    wo.endQuerySequence();
    wo.nextSequence(1, 0);
    wo.nextQuery(BidirectionalFrame.FORWARD, 1);
    wo.setPosition(0);
    wo.endPosition();
    wo.endQuery();
    wo.endQuerySequence();
    wo.endAll();
    assertEquals(HEADER + EXPECTED2, out.toString());
    assertEquals(4.0, wo.score());
    assertEquals("1" + StringUtils.LS, unmappedOut.toString());
  }

  public final void testToStringStringBuilder() {
    final Appendable out = new StringWriter();
    final Appendable unmappedOut = new StringWriter();
    assertEquals("WordOutput", new WordOutput(UnidirectionalFrame.values(), out, unmappedOut).toString());
  }

}
