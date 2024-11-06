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
