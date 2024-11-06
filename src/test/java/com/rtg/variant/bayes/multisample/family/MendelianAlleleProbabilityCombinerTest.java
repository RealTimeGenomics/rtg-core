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
    assertEquals(Math.log(0.00001 / 6), c.probabilityLn(code, code.code(0, 1), code.code(0, 0), code.code(1, 1)), 1e-8);
  }
  public void testDenovoProbabilityLotsOfAlleles() {
    final MendelianAlleleProbabilityCombiner c = new MendelianAlleleProbabilityCombiner(MendelianAlleleProbabilityDiploid.SINGLETON, MendelianAlleleProbabilityDiploidDeNovo.SINGLETON, Math.log(0.00002), Math.log(0.00001), 0);
    final CodeDiploid code = new CodeDiploid(4);
    assertEquals(Math.log(0.00001 / 12), c.probabilityLn(code, code.code(0, 1), code.code(1, 2), code.code(2, 3)), 1e-8);
  }

  public void testDenovo() {
    final MendelianAlleleProbabilityCombiner c = new MendelianAlleleProbabilityCombiner(MendelianAlleleProbabilityDiploid.SINGLETON, MendelianAlleleProbabilityDiploidDeNovo.SINGLETON, Math.log(0.0000025), Math.log(0.00001), 0);
    final CodeDiploid code = new CodeDiploid(4);
    assertTrue(c.isDenovo(code, code.code(0, 0), code.code(0, 0), code.code(0, 1)));
    assertTrue(c.isDenovo(code, code.code(1, 1), code.code(1, 1), code.code(0, 1)));
  }

}
