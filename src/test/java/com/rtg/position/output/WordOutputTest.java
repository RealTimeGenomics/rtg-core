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
