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

package com.rtg.calibrate;

import com.rtg.sam.BadSuperCigarException;

import htsjdk.samtools.SAMRecord;

import junit.framework.TestCase;

/**
 */
public class CovariateBaseQualityTest extends TestCase {


  public void test() throws BadSuperCigarException {
    final Covariate var = new CovariateBaseQuality();
    assertEquals(64, var.size());
    assertEquals(64, var.newSize());
    final SAMRecord sam = new SAMRecord(null);
    sam.setBaseQualityString("ABCDEFGHIJ");

    final Calibrator calwtf = new Calibrator(new Covariate[] {new CovariateReadGroup()}, null);
    final CalibratorCigarParser parser = new CalibratorCigarParser(calwtf) {
      @Override
      public int getCurrentQuality() {
        return 'J' - '!';
      }
    };
    assertEquals('J' - '!', var.value(sam, parser));
    assertEquals("14", var.valueString(14));
    assertEquals(14, var.parse("14"));
  }

  public void testNoQuality() throws BadSuperCigarException {
    final Covariate var = new CovariateBaseQuality();
    final SAMRecord sam = new SAMRecord(null);
    sam.setBaseQualityString("*");

    final Calibrator calwtf = new Calibrator(new Covariate[] {new CovariateReadGroup()}, null);
    final CalibratorCigarParser parser = new CalibratorCigarParser(calwtf) {
      @Override
      public int getReadPosition() {
        return 0;
      }
    };
    assertEquals('5' - '!', var.value(sam, parser));
  }

  public void testNameType() {
    final Covariate var = new CovariateBaseQuality();
    assertEquals("basequality", var.name());
    assertEquals(CovariateEnum.BASEQUALITY, var.getType());
  }
}
