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

import static com.rtg.util.StringUtils.LS;
import static com.rtg.util.StringUtils.TAB;

import java.util.Arrays;

import com.rtg.variant.match.AlignmentMatch;

import junit.framework.TestCase;

/**
 */
public class MatchMultiSetTest extends TestCase {

  public void test() {
    final MatchMultiSet mms = new MatchMultiSet();
    mms.integrity();
    final String exp = ""
      + "[ MatchMultiSet size=0" + LS
      + "]" + LS
      ;
    assertEquals(exp, mms.toString());
    assertEquals(0, mms.size());
    assertEquals(0, mms.totalCount());
    assertEquals("[]", Arrays.toString(mms.names()));
    final String out = mms.output(false);
    assertEquals("", out);
  }

  public void test1() {
    final MatchMultiSet mms = new MatchMultiSet("");
    mms.integrity();
    final String exp = ""
      + "[ MatchMultiSet size=1" + LS
      + " > 0:0.000" + LS
      + "]" + LS
      ;
    assertEquals(exp, mms.toString());
    assertEquals(1, mms.size());
    assertEquals(0, mms.totalCount());
    assertEquals(1, mms.names().length);
    assertEquals("[]", Arrays.toString(mms.names()));
    final String out = mms.output(false);
    assertEquals("", out);
  }

  public void test2() {
    final MatchMultiSet mms = new MatchMultiSet("");
    mms.integrity();
    mms.add(new AlignmentMatch(null, "", null, 20, 0, 0, 20), 0.01);
    mms.add(new AlignmentMatch(null, "ACTG", null, 20, 0, 4, 20), 0.02);
    final String exp = ""
      + "[ MatchMultiSet size=2" + LS
      + " > 1:0.010" + LS
      + "ACTG > 1:0.020" + LS
      + "]" + LS
      ;
    assertEquals(exp, mms.toString());
    assertEquals(2, mms.size());
    assertEquals(2, mms.totalCount());
    assertEquals("[, ACTG]", Arrays.toString(mms.names()));
    final String out = mms.output(false);
    assertEquals(TAB + "" + TAB + "1" + TAB + "0.010" + TAB + "ACTG" + TAB + "1" + TAB + "0.020", out);
    mms.add(new AlignmentMatch(null, "ACTG", null, 20, 0, 4, 20), 0.02);
    assertEquals(2, mms.size());
    assertEquals(3, mms.totalCount());
    mms.reset();

    assertEquals(0, mms.size());
    assertEquals(0, mms.totalCount());
  }

  public void test2lc() {
    final MatchMultiSet mms = new MatchMultiSet("");
    mms.integrity();
    mms.add(new AlignmentMatch(null, "", null, 20, 0, 0, 20), 0.01);
    mms.add(new AlignmentMatch(null, "ACTG", null, 20, 0, 4, 20), 0.02);
    final String exp = ""
      + "[ MatchMultiSet size=2" + LS
      + " > 1:0.010" + LS
      + "ACTG > 1:0.020" + LS
      + "]" + LS
      ;
    assertEquals(exp, mms.toString());
    assertEquals(2, mms.size());
    assertEquals(2, mms.totalCount());
    assertEquals("[, ACTG]", Arrays.toString(mms.names()));
    final String out = mms.output(true);
    assertEquals(TAB + "" + TAB + "1" + TAB + "0.010" + TAB + "actg" + TAB + "1" + TAB + "0.020", out);
    mms.add(new AlignmentMatch(null, "ACTG", null, 20, 0, 4, 20), 0.02);
    assertEquals(2, mms.size());
    assertEquals(3, mms.totalCount());
  }
}
