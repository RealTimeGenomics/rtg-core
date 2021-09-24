/*
 * Copyright (c) 2017. Real Time Genomics Limited.
 *
 * Use of this source code is bound by the Real Time Genomics Limited Software Licence Agreement
 * for Academic Non-commercial Research Purposes only.
 *
 * If you did not receive a license accompanying this file, a copy must first be obtained by email
 * from support@realtimegenomics.com.  On downloading, using and/or continuing to use this source
 * code you accept the terms of that license agreement and any amendments to those terms that may
 * be made from time to time by Real Time Genomics Limited.
 */
package com.rtg.variant.bayes.complex;

import static com.rtg.variant.bayes.complex.AbstractDecomposerTest.createVariant;

import com.rtg.variant.Variant;

import junit.framework.TestCase;

/**
 */
public class TrimmerTest extends TestCase {

  Decomposer getDecomposer() {
    return new Trimmer(SimpleDecomposerTest.TRIGGER);
  }

  public void testTrimSimpleSnps() {
    //homozygous
    //                               0123456789
    final Variant v = createVariant("ACGTCTGTCT", "ACCTCTGTCT", null);
    final Variant trim = getDecomposer().decompose(v).get(0);
    assertEquals('C', trim.getLocus().getPreviousRefNt());
    assertEquals("G", trim.getLocus().getRefNts());
    assertEquals(2, trim.getLocus().getStart());
    assertEquals(3, trim.getLocus().getEnd());
    assertEquals("C", trim.getSample(0).getName());

    //heterozygous
    //                                0123456789
    final Variant v2 = createVariant("ACGTCTGTCT", "ACCTCTGTCT", "ACGTCTGTCT");
    final Variant trim2 = getDecomposer().decompose(v2).get(0);
    assertEquals('C', trim2.getLocus().getPreviousRefNt());
    assertEquals("G", trim2.getLocus().getRefNts());
    assertEquals(2, trim2.getLocus().getStart());
    assertEquals(3, trim2.getLocus().getEnd());
    assertEquals("C:G", trim2.getSample(0).getName());
  }

  public void testTrimDelete() {
    //homozygous delete
    //                                               0123456789
    final Variant vHomozygousDelete = createVariant("ACGTCTGTCT", "ACGTCTGT", null);
    final Variant trim = getDecomposer().decompose(vHomozygousDelete).get(0);
    assertEquals('G', trim.getLocus().getPreviousRefNt());
    assertEquals("TC", trim.getLocus().getRefNts());
    assertEquals(7, trim.getLocus().getStart());
    assertEquals(9, trim.getLocus().getEnd());
    assertEquals("", trim.getSample(0).getName());

    //heterozygous deletes
    final Variant vHeterozygousDelete = createVariant("ACGTCTGTCT", "ACGTCTGT", "ACGTCTG");
    final Variant trim2 = getDecomposer().decompose(vHeterozygousDelete).get(0);
    assertEquals('G', trim2.getLocus().getPreviousRefNt());
    assertEquals("TCT", trim2.getLocus().getRefNts());
    assertEquals(7, trim2.getLocus().getStart());
    assertEquals(10, trim2.getLocus().getEnd());
    assertEquals("T:", trim2.getSample(0).getName());

  //heterozygous deletes
    final Variant vHeterozygousDelete2 = createVariant("ACGTCTGTCT", "ACGTCTG", "ACGTCTGT");
    final Variant trim3 = getDecomposer().decompose(vHeterozygousDelete2).get(0);
    assertEquals('G', trim3.getLocus().getPreviousRefNt());
    assertEquals("TCT", trim3.getLocus().getRefNts());
    assertEquals(":T", trim3.getSample(0).getName());
  }

  public void testTrimInsert() {
    //homozygous
    //                                               0123456789
    final Variant vHomozygousInsert = createVariant("ACGTCTGTCT", "ACGTTTCTGTCT", null);
    final Variant trim = getDecomposer().decompose(vHomozygousInsert).get(0);
    assertEquals('G', trim.getLocus().getPreviousRefNt());
    assertEquals("", trim.getLocus().getRefNts());
    assertEquals(3, trim.getLocus().getStart());
    assertEquals(3, trim.getLocus().getEnd());
    assertEquals("TT", trim.getSample(0).getName());

    //heterozygous
    //                                               0123456789
    final Variant vHeterozyygousInsert = createVariant("ACGTCTGTCT", "ACGTTTCTGTCT", "ACGTCCCTGTCT");
    final Variant trim2 = getDecomposer().decompose(vHeterozyygousInsert).get(0);
    assertEquals('T', trim2.getLocus().getPreviousRefNt());
    assertEquals("", trim2.getLocus().getRefNts());
    assertEquals(4, trim2.getLocus().getStart());
    assertEquals(4, trim2.getLocus().getEnd());
    assertEquals("TT:CC", trim2.getSample(0).getName());
  }

  public void testTrimComplex() {
  //homozygous
    //                                               0123456789
    final Variant vHomozygousIndel = createVariant("ACGTCTGTCT", "ACGTTTCTGGCT", null);
    final Variant trim = getDecomposer().decompose(vHomozygousIndel).get(0);
    assertEquals('T', trim.getLocus().getPreviousRefNt());
    assertEquals("CTGT", trim.getLocus().getRefNts());
    assertEquals(4, trim.getLocus().getStart());
    assertEquals(8, trim.getLocus().getEnd());
    assertEquals("TTCTGG", trim.getSample(0).getName());

    //heterozygous
    // 01|234567---8|9
    // AC|GTCTAT---C|T
    // AC|--CTATTTTC|T
    // AC|GTCTA----G|T
    //                                                  0123456789
    final Variant vHeterozyygousIndel = createVariant("ACGTCTATCT", "ACCTATTTTCT", "ACGTCTAGT");
    final Variant trim2 = getDecomposer().decompose(vHeterozyygousIndel).get(0);
    assertEquals('C', trim2.getLocus().getPreviousRefNt());
    assertEquals("GTCTATC", trim2.getLocus().getRefNts());
    assertEquals(2, trim2.getLocus().getStart());
    assertEquals(9, trim2.getLocus().getEnd());
    assertEquals("CTATTTTC:GTCTAG", trim2.getSample(0).getName());
  }

  public void testTrimRepeat() {
    //homozygous
    //
    final Variant vHomozygousIndel = createVariant("TA", "TAA", null);
    final Variant trim = getDecomposer().decompose(vHomozygousIndel).get(0);
    assertEquals('T', trim.getLocus().getPreviousRefNt());
    assertEquals("", trim.getLocus().getRefNts());
    assertEquals(1, trim.getLocus().getStart());
    assertEquals(1, trim.getLocus().getEnd());
    assertEquals("A", trim.getSample(0).getName());

    //heterozygous
    // AC|GTCTA----G|T
    final Variant vHeterozyygousIndel = createVariant("TA", "TAA", "TA");
    final Variant trim2 = getDecomposer().decompose(vHeterozyygousIndel).get(0);
    assertEquals('T', trim2.getLocus().getPreviousRefNt());
    assertEquals("", trim2.getLocus().getRefNts());
    assertEquals(1, trim2.getLocus().getStart());
    assertEquals(1, trim2.getLocus().getEnd());
    assertEquals("A:", trim2.getSample(0).getName());
  }


}
