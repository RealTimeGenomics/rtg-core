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

import java.util.List;

import com.rtg.variant.Variant;
import com.rtg.variant.VariantLocus;
import com.rtg.variant.bayes.multisample.VariantAlleleTrigger;

/**
 */
public class AligningDecomposerTest extends AbstractDecomposerTest {

  @Override
  public Decomposer getDecomposer() {
    return new AligningDecomposer();
  }

  @Override
  public Decomposer getDecomposer(final DenovoChecker denovoChecker) {
    return new AligningDecomposer(denovoChecker, SimpleDecomposerTest.TRIGGER);
  }

  @Override
  public Decomposer getDecomposer(final VariantAlleleTrigger variantAlleleTrigger) {
    return new AligningDecomposer(null, variantAlleleTrigger);
  }

  public void testReference() {
    // homozygous identity
    final Variant v = createVariant("GTCCTAACT", "GTCCTAACT", null);
    assertTrue(getDecomposer().decompose(v).isEmpty());

    // heterozygous identity
    final Variant v2 = createVariant("GTCCTAACT", "GTCCTAACT", "GTCCTAACT");
    assertTrue(getDecomposer().decompose(v2).isEmpty());
  }

  public void testHomozygousDifferentLengths() {
    final List<Variant> list = getDecomposer().decompose(createVariant("GTCCTAACT", "GTCCAACT", null));
    assertEquals(1, list.size());
    final Variant first = list.get(0);
    assertEquals("", first.getSample(0).getName());
    assertEquals(30.5, first.getSample(0).getPosterior());
    assertEquals(false, first.getSample(0).isIdentity());
    final VariantLocus locus = first.getLocus();
    assertEquals("T", locus.getRefNts());
    assertEquals("chr", locus.getSequenceName());
    assertEquals(4, locus.getStart());
    assertEquals(5, locus.getEnd());
    assertEquals('C', locus.getPreviousRefNt());
  }

  public void testHeterozygousDifferentLengths() {
    final List<Variant> list = getDecomposer().decompose(createVariant("GTCCTAACT", "GTCCTAA", "GTCCTAACTAA"));
    assertEquals(1, list.size());
    final Variant first = list.get(0);
    assertEquals(":CTAA", first.getSample(0).getName());
    assertEquals(30.5, first.getSample(0).getPosterior());
    assertEquals(false, first.getSample(0).isIdentity());
    final VariantLocus locus = first.getLocus();
    assertEquals("CT", locus.getRefNts());
    assertEquals("chr", locus.getSequenceName());
    assertEquals(7, locus.getStart());
    assertEquals(9, locus.getEnd());
    assertEquals('A', locus.getPreviousRefNt());
  }

  public void testTrimIdentity() {
    assertTrue(getDecomposer().decompose(createVariant("ACC", "ACC", null)).isEmpty());
    assertTrue(getDecomposer().decompose(createVariant("ACC", "ACC", "ACC")).isEmpty());
  }

  @Override
  public void testTrimInsertHetero() {
    final List<Variant> decomposed = getDecomposer().decompose(createVariant("ACGTCTGTCT", "ACGTTTCTGTCT", "ACGTCCCTGTCT"));
    final Variant v0 = decomposed.get(0);
    assertEquals("TT:", v0.getSample(0).getName());
    assertEquals(3, v0.getLocus().getStart());
    assertEquals(3, v0.getLocus().getEnd());
    assertEquals('G', v0.getLocus().getPreviousRefNt());
    assertEquals("", v0.getLocus().getRefNts());
    final Variant v1 = decomposed.get(1);
    assertEquals(":CC", v1.getSample(0).getName());
    assertEquals(4, v1.getLocus().getStart());
    assertEquals(4, v1.getLocus().getEnd());
    assertEquals('T', v1.getLocus().getPreviousRefNt());
    assertEquals("", v1.getLocus().getRefNts());
  }

  public void testTrimComplexHomo() {
    final List<Variant> decompose = getDecomposer().decompose(createVariant("ACGTCTGTCT", "ACGTTTCTGGCT", null));
    assertEquals(2, decompose.size());
    final Variant v0 = decompose.get(0);
    assertEquals("TT", v0.getSample(0).getName());
    assertEquals("", v0.getLocus().getRefNts());
    assertEquals('G', v0.getLocus().getPreviousRefNt());
    assertEquals(3, v0.getLocus().getStart());
    assertEquals(3, v0.getLocus().getEnd());
    final Variant v1 = decompose.get(1);
    assertEquals("G", v1.getSample(0).getName());
    assertEquals("T", v1.getLocus().getRefNts());
    assertEquals('G', v1.getLocus().getPreviousRefNt());
    assertEquals(7, v1.getLocus().getStart());
    assertEquals(8, v1.getLocus().getEnd());
  }

  public void testTrimComplexHetero() {
    // 01 2 3 4  5 678 9
    // AC G T C  T ATC T
    // AC C T AT T TTC T
    // AC G T C  T AG  T
    final List<Variant> decompose = getDecomposer().decompose(createVariant("ACGTCTATCT", "ACCTATTTTCT", "ACGTCTAGT"));
    assertEquals(3, decompose.size());
    final Variant v0 = decompose.get(0);
    assertEquals("C:G", v0.getSample(0).getName());
    assertEquals("G", v0.getLocus().getRefNts());
    assertEquals('C', v0.getLocus().getPreviousRefNt());
    assertEquals(2, v0.getLocus().getStart());
    assertEquals(3, v0.getLocus().getEnd());
    final Variant v1 = decompose.get(1);
    assertEquals("AT:C", v1.getSample(0).getName());
    assertEquals("C", v1.getLocus().getRefNts());
    assertEquals('T', v1.getLocus().getPreviousRefNt());
    assertEquals(4, v1.getLocus().getStart());
    assertEquals(5, v1.getLocus().getEnd());
    final Variant v2 = decompose.get(2);
    assertEquals("TTC:AG", v2.getSample(0).getName());
    assertEquals("ATC", v2.getLocus().getRefNts());
    assertEquals('T', v2.getLocus().getPreviousRefNt());
    assertEquals(6, v2.getLocus().getStart());
    assertEquals(9, v2.getLocus().getEnd());
  }
}
