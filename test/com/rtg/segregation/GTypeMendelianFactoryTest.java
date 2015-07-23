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
import com.rtg.segregation.GTypeMendelianFactory.Swapped;
import com.rtg.segregation.GTypeMendelianFactory.TriplePloidy;

import junit.framework.TestCase;

/**
 */
public class GTypeMendelianFactoryTest extends TestCase {

  private static GType sA;
  private static GType sB;
  private static GType sC;
  private static GType sD;
  private static GType sE;
  private static GType sF;
  private static GType sG;
  private static GType sH;
  private static GType sI;
  private static GType sJ;
  private static GType sK;
  static {
    try {
      sA = new GType("0/0", Ploidy.DIPLOID);
      sB = new GType("0/1", Ploidy.DIPLOID);
      sC = new GType("1/1", Ploidy.DIPLOID);
      sD = new GType("0/2", Ploidy.DIPLOID);
      sE = new GType("1/2", Ploidy.DIPLOID);
      sF = new GType("2/2", Ploidy.DIPLOID);
      sG = new GType("0", Ploidy.HAPLOID);
      sH = new GType("1", Ploidy.HAPLOID);
      sI = new GType(".", Ploidy.NONE);
      sJ = new GType("0", Ploidy.POLYPLOID);
      sK = new GType("1", Ploidy.POLYPLOID);
    } catch (MismatchingPloidyException e) { }
  }

  public void testDDD() {
    final GTypeMendelian checker = GTypeMendelianFactory.getGTypeMendelian(Ploidy.DIPLOID, Ploidy.DIPLOID, Ploidy.DIPLOID);
    assertNotNull(checker);
    assertTrue(checker.isMendelian(sA, sA, sA));
    assertFalse(checker.isMendelian(sA, sA, sB));
    assertFalse(checker.isMendelian(sA, sC, sF));
    assertFalse(checker.isMendelian(sA, sB, sC));
    assertFalse(checker.isMendelian(sD, sA, sE));
    assertFalse(checker.isMendelian(sE, sC, sB));
    assertTrue(checker.isMendelian(sC, sA, sB));
    assertTrue(checker.isMendelian(sC, sB, sC));
    assertTrue(checker.isMendelian(sB, sE, sE));
    assertTrue(checker.isMendelian(sD, sB, sE));
  }

  public void testHDD() {
    final GTypeMendelian checker = GTypeMendelianFactory.getGTypeMendelian(Ploidy.HAPLOID, Ploidy.DIPLOID, Ploidy.DIPLOID);
    assertNotNull(checker);
    assertTrue(checker.isMendelian(sG, sA, sA));
    assertFalse(checker.isMendelian(sH, sA, sA));
    assertFalse(checker.isMendelian(sG, sC, sA));
    assertTrue(checker.isMendelian(sH, sA, sB));
  }

  public void testHDH() {
    final GTypeMendelian checker = GTypeMendelianFactory.getGTypeMendelian(Ploidy.HAPLOID, Ploidy.DIPLOID, Ploidy.HAPLOID);
    assertNotNull(checker);
    assertTrue(checker.isMendelian(sG, sA, sG));
    assertFalse(checker.isMendelian(sG, sA, sH));
  }

  public void testNHN() {
    final GTypeMendelian checker = GTypeMendelianFactory.getGTypeMendelian(Ploidy.NONE, Ploidy.HAPLOID, Ploidy.NONE);
    assertNotNull(checker);
    assertTrue(checker.isMendelian(sI, sG, sI));
  }

  public void testHNH() {
    final GTypeMendelian checker = GTypeMendelianFactory.getGTypeMendelian(Ploidy.HAPLOID, Ploidy.NONE, Ploidy.HAPLOID);
    assertNotNull(checker);
    assertTrue(checker.isMendelian(sG, sI, sG));
    assertFalse(checker.isMendelian(sG, sI, sH));
  }

  public void testPPP() {
    final GTypeMendelian checker = GTypeMendelianFactory.getGTypeMendelian(Ploidy.POLYPLOID, Ploidy.POLYPLOID, Ploidy.POLYPLOID);
    assertNotNull(checker);
    assertTrue(checker.isMendelian(sJ, sJ, sJ));
    assertFalse(checker.isMendelian(sK, sJ, sK));
  }

  public void testSwapped() {
    final GTypeMendelian swp = new Swapped(new GTypeMendelian() {
      @Override
      public boolean isMendelian(GType fa, GType mo, GType ch) {
        return fa.a() == ch.a();
      }
    });
    assertTrue(swp.isMendelian(sC, sA, sA));
    assertFalse(swp.isMendelian(sA, sC, sA));
  }

  public void testTriplePloidy() {
    final TriplePloidy base = new TriplePloidy(Ploidy.DIPLOID, Ploidy.HAPLOID, Ploidy.POLYPLOID);
    final TriplePloidy one = new TriplePloidy(Ploidy.HAPLOID, Ploidy.HAPLOID, Ploidy.POLYPLOID);
    final TriplePloidy two = new TriplePloidy(Ploidy.DIPLOID, Ploidy.DIPLOID, Ploidy.POLYPLOID);
    final TriplePloidy three = new TriplePloidy(Ploidy.DIPLOID, Ploidy.HAPLOID, Ploidy.HAPLOID);
    assertFalse(base.equals(null));
    assertTrue(base.equals(base));
    assertFalse(base.equals(one));
    assertFalse(base.equals(two));
    assertFalse(base.equals(three));
    assertFalse(base.hashCode() == one.hashCode());
  }
}
