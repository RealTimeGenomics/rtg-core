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

package com.rtg.variant;

import junit.framework.TestCase;

/**
 */
public class VariantSampleInfoTest extends TestCase {

  public void testHowItWorks() {
    final VariantSampleInfo info = new VariantSampleInfo("Blah");
    assertEquals("Blah", info.stringValue());
    assertNull(info.doubleValue());
    assertNull(info.intValue());
    assertFalse(info.boolValue());  //note this behaviour!

    final VariantSampleInfo info2 = new VariantSampleInfo(1.5);
    assertEquals(1.5, info2.doubleValue());
    assertNull(info2.stringValue());

    final VariantSampleInfo info3 = new VariantSampleInfo(true);
    assertTrue(info3.boolValue());

  }

  public void testEnum() {
    com.rtg.util.TestUtils.testEnum(VariantSampleInfo.VariantFormatEnum.class, "[COVERAGE, COVERAGE_CORRECTION, AMBIGUITY_RATIO, NONIDENTITY_POSTERIOR, STATISTICS, HOEFFDING_STRAND_BIAS_A1, HOEFFDING_STRAND_BIAS_A2, HOEFFDING_ALLELE_BALANCE_HET, HOEFFDING_ALLELE_BALANCE_HOM, HOEFFDING_READ_POSITION, HOEFFDING_UNMATED_BIAS_A1, HOEFFDING_UNMATED_BIAS_A2, UNPLACED_RATIO, SOMATIC_SCORE]");
  }

}
