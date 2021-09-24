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

import java.util.ArrayList;
import java.util.List;

import com.rtg.variant.match.Match;

import junit.framework.TestCase;

/**
 */
public class DescriptionComplexTest extends TestCase {

  public void testEmpty() {
    final List<Match> ml = new ArrayList<>();
    final DescriptionComplex cmpx = new DescriptionComplex(ml);
    assertEquals(0, cmpx.size());
  }

  public void test() {
    final List<Match> ml = new ArrayList<>();
    ml.add(HypothesesComplexTest.match("A", 20));
    final Match matchAA = HypothesesComplexTest.match("AA", 20);
    ml.add(matchAA);
    ml.add(HypothesesComplexTest.match("AAA", 20));
    ml.add(HypothesesComplexTest.match("AAAA", 20));
    ml.add(HypothesesComplexTest.match("AAAAA", 20));
    final DescriptionComplex cmpx = new DescriptionComplex(ml);
    assertEquals(5, cmpx.size());
    assertEquals("AA", cmpx.name(1));
    assertEquals(matchAA, cmpx.match(1));
    assertEquals(1, cmpx.minLength());
    assertEquals(5, cmpx.maxLength());
  }

  public void test0() {
    final List<Match> ml = new ArrayList<>();
    final DescriptionComplex cmpx = new DescriptionComplex(ml);
    assertEquals(0, cmpx.size());
    assertEquals(Integer.MAX_VALUE, cmpx.minLength());
    assertEquals(0, cmpx.maxLength());
  }

}
