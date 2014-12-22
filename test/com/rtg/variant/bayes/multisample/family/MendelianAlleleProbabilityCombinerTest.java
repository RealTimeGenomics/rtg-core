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

}
