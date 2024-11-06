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
