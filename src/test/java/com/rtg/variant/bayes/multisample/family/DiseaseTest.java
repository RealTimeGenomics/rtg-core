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
package com.rtg.variant.bayes.multisample.family;

import java.io.BufferedReader;
import java.io.StringReader;

import com.rtg.relation.Family;
import com.rtg.relation.FamilyTest;
import com.rtg.relation.RelationshipsFileParser;
import com.rtg.variant.GenomePriorParams;
import com.rtg.variant.bayes.Code;
import com.rtg.variant.bayes.Hypotheses;
import com.rtg.variant.bayes.snp.HypothesesSnp;
import com.rtg.variant.util.arithmetic.SimplePossibility;

import junit.framework.TestCase;

/**
 */
public class DiseaseTest extends TestCase {

  public void testClassic() throws Exception {
    final Family f = Family.getFamily(RelationshipsFileParser.load(new BufferedReader(new StringReader(FamilyTest.RELATIONS))));
    assertTrue(f.isDiseased("childb"));
    assertTrue(f.isDiseased(f.getMother()));
    final Hypotheses<?> hypotheses = new HypothesesSnp(SimplePossibility.SINGLETON, GenomePriorParams.builder().create(), false, 0);
    final Code code = hypotheses.code();
    final Disease d = new Disease(f, hypotheses, code.code(0), code.code(0, 3), code.code(0), code.code(3));
    assertEquals(4, d.explanation(0));
    assertEquals(4, d.explanation(1));
    assertEquals(4, d.explanation(2));
    assertEquals(0, d.explanation(3));
  }

  public void testClassic2() throws Exception {
    final Family f = Family.getFamily(RelationshipsFileParser.load(new BufferedReader(new StringReader(FamilyTest.RELATIONS))));
    assertTrue(f.isDiseased("childb"));
    assertTrue(f.isDiseased(f.getMother()));
    final Hypotheses<?> hypotheses = new HypothesesSnp(SimplePossibility.SINGLETON, GenomePriorParams.builder().create(), false, 0);
    final Code code = hypotheses.code();
    final Disease d = new Disease(f, hypotheses, code.code(3), code.code(3, 0), code.code(3), code.code(0));
    assertEquals(0, d.explanation(0));
    assertEquals(1, d.explanation(1));
    assertEquals(1, d.explanation(2));
    assertEquals(1, d.explanation(3));
  }

  public void testNoExplanation() throws Exception {
    final Family f = Family.getFamily(RelationshipsFileParser.load(new BufferedReader(new StringReader(FamilyTest.RELATIONS))));
    final Hypotheses<?> hypotheses = new HypothesesSnp(SimplePossibility.SINGLETON, GenomePriorParams.builder().create(), false, 0);
    final Code code = hypotheses.code();
    final Disease d = new Disease(f, hypotheses, code.code(1, 2), code.code(2, 3), code.code(1, 1), code.code(1, 2));
    assertEquals(0, d.explanation(0));
    assertEquals(0, d.explanation(1));
    assertEquals(0, d.explanation(2));
    assertEquals(0, d.explanation(3));
  }

}
