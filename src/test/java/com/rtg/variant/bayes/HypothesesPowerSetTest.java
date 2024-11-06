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

package com.rtg.variant.bayes;

import com.rtg.reference.Ploidy;
import com.rtg.variant.bayes.snp.DescriptionCommon;
import com.rtg.variant.util.arithmetic.LogPossibility;

import junit.framework.TestCase;

/**
 * Tests the corresponding class.
 */
public class HypothesesPowerSetTest extends TestCase {

  public void test() {
    final DescriptionCommon desc = new DescriptionCommon("A", "C", "G", "T");
    final HypothesesPowerSet<DescriptionCommon> hyp = new HypothesesPowerSet<>(desc, LogPossibility.SINGLETON, 2);
    assertEquals(15, hyp.size());
    assertFalse(hyp.valid(-1));
    assertTrue(hyp.valid(0));
    assertTrue(hyp.valid(14));
    assertFalse(hyp.valid(15));
    assertEquals("A:C:G:T", hyp.name(14));
    assertEquals("G", hyp.name(hyp.reference()));
    assertEquals(LogPossibility.SINGLETON, hyp.arithmetic());
    assertEquals(7, hyp.maxNameLength());
    assertEquals(Ploidy.POLYPLOID, hyp.ploidy());
    assertTrue(hyp.code() instanceof CodePowerSet);
  }
}
