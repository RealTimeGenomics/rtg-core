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

package com.rtg.segregation;

import com.rtg.reference.Ploidy;

import junit.framework.TestCase;

/**
 */
public class GTypeTest extends TestCase {

  public void test() throws MismatchingPloidyException {
    final GType gt = new GType("0/1", Ploidy.DIPLOID);
    gt.integrity();
    assertEquals("0_1", gt.toString());
    assertEquals(0, gt.a());
    assertEquals(1, gt.b());
    assertFalse(gt.isSingleAllele());
  }

  public void testIsHomozygous() throws MismatchingPloidyException {
    final GType gt = new GType("0|0", Ploidy.DIPLOID);
    gt.integrity();
    assertEquals("0_0", gt.toString());
    assertEquals(0, gt.a());
    assertEquals(0, gt.b());
    assertTrue(gt.isSingleAllele());
  }

  public void testIsMendelian() throws MismatchingPloidyException {
    check("0/0", "0|0", "0/0", true);
    check("0|0", "0|1", "0|0", true);
    check("0/0", "0/1", "0/1", true);
    check("0/0", "1/1", "0/1", true);

    check("0/1", "2/3", "0/2", true);
    check("0/1", "2/3", "0/3", true);
    check("0/1", "2/3", "1/2", true);
    check("0/1", "2/3", "1/3", true);
    check("0/2", "1/3", "1/2", true);


    check("0/0", "0/0", "0/1", false);

  }

  private void check(String fa, String mo, String ch, boolean exp) throws MismatchingPloidyException {
    final GType gfa = new GType(fa, Ploidy.DIPLOID);
    final GType gmo = new GType(mo, Ploidy.DIPLOID);
    final GType gch = new GType(ch, Ploidy.DIPLOID);
    final GTypeMendelian mendelianChecker = GTypeMendelianFactory.getGTypeMendelian(Ploidy.DIPLOID, Ploidy.DIPLOID, Ploidy.DIPLOID);
    assertEquals(exp, mendelianChecker.isMendelian(gfa, gmo, gch));
    assertEquals(exp, mendelianChecker.isMendelian(gmo, gfa, gch));
  }
}
