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
