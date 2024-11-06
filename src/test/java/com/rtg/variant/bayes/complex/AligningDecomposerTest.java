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
