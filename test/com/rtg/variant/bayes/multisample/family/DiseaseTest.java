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
