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

import com.rtg.variant.bayes.CodeDiploid;

import junit.framework.TestCase;

/**
 */
public class MendelianAlleleProbabilityCombinerTest extends TestCase {

  public void test() {
    final MendelianAlleleProbabilityCombiner c = new MendelianAlleleProbabilityCombiner(MendelianAlleleProbabilityDiploid.SINGLETON, MendelianAlleleProbabilityDiploidDeNovo.SINGLETON, Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY, 0);
    final CodeDiploid code = new CodeDiploid(4);
    assertEquals(MendelianAlleleProbabilityDiploid.SINGLETON.probabilityLn(code, 0, 1, 1), c.probabilityLn(code, 0, 1, 1));
  }

  public void testDenovoProbabilityRef() {
    final MendelianAlleleProbabilityCombiner c = new MendelianAlleleProbabilityCombiner(MendelianAlleleProbabilityDiploid.SINGLETON, MendelianAlleleProbabilityDiploidDeNovo.SINGLETON, Math.log(0.00002), Math.log(0.00001), 0);
    final CodeDiploid code = new CodeDiploid(4);
    assertEquals(Math.log(0.00002 / 6), c.probabilityLn(code, code.code(0, 0), code.code(0, 0), code.code(0, 1)), 1e-8);
  }
  public void testDenovoProbabilityNonRef() {
    final MendelianAlleleProbabilityCombiner c = new MendelianAlleleProbabilityCombiner(MendelianAlleleProbabilityDiploid.SINGLETON, MendelianAlleleProbabilityDiploidDeNovo.SINGLETON, Math.log(0.00002), Math.log(0.00001), 0);
    final CodeDiploid code = new CodeDiploid(4);
    assertEquals(Math.log(0.00001 / 6), c.probabilityLn(code, code.code(1, 1), code.code(1, 1), code.code(1, 2)), 1e-8);
  }

  public void testDenovoProbabilityNonRefDiploidParent() {
    final MendelianAlleleProbabilityCombiner c = new MendelianAlleleProbabilityCombiner(MendelianAlleleProbabilityDiploid.SINGLETON, MendelianAlleleProbabilityDiploidDeNovo.SINGLETON, Math.log(0.00002), Math.log(0.00001), 0);
    final CodeDiploid code = new CodeDiploid(4);
    assertEquals(Math.log(0.00001/ 6), c.probabilityLn(code, code.code(0, 1), code.code(0, 0), code.code(1, 1)), 1e-8);
  }
  public void testDenovoProbabilityLotsOfAlleles() {
    final MendelianAlleleProbabilityCombiner c = new MendelianAlleleProbabilityCombiner(MendelianAlleleProbabilityDiploid.SINGLETON, MendelianAlleleProbabilityDiploidDeNovo.SINGLETON, Math.log(0.00002), Math.log(0.00001), 0);
    final CodeDiploid code = new CodeDiploid(4);
    assertEquals(Math.log(0.00001/ 12), c.probabilityLn(code, code.code(0, 1), code.code(1, 2), code.code(2, 3)), 1e-8);
  }

  public void testDenovo() {
    final MendelianAlleleProbabilityCombiner c = new MendelianAlleleProbabilityCombiner(MendelianAlleleProbabilityDiploid.SINGLETON, MendelianAlleleProbabilityDiploidDeNovo.SINGLETON, Math.log(0.0000025), Math.log(0.00001), 0);
    final CodeDiploid code = new CodeDiploid(4);
    assertTrue(c.isDenovo(code, code.code(0, 0), code.code(0, 0), code.code(0, 1)));
    assertTrue(c.isDenovo(code, code.code(1, 1), code.code(1,1), code.code(0, 1)));
  }

}
