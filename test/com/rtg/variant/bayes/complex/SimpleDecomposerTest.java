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
package com.rtg.variant.bayes.complex;

import java.util.List;

import com.rtg.reference.Ploidy;
import com.rtg.variant.Variant;
import com.rtg.variant.VariantLocus;
import com.rtg.variant.VariantSample;
import com.rtg.variant.bayes.multisample.VariantAlleleTrigger;
import com.rtg.variant.bayes.snp.DescriptionCommon;
import com.rtg.variant.format.VariantOutputVcfFormatterTest;

/**
 */
public class SimpleDecomposerTest extends AbstractDecomposerTest {

  static final VariantAlleleTrigger TRIGGER = new VariantAlleleTrigger(0, 0.0);

  @Override
  public Decomposer getDecomposer() {
    return new SimpleDecomposer();
  }

  @Override
  public Decomposer getDecomposer(final DenovoChecker denovoChecker) {
    return new SimpleDecomposer(denovoChecker, TRIGGER);
  }

  @Override
  public Decomposer getDecomposer(final VariantAlleleTrigger variantAlleleTrigger) {
    return new SimpleDecomposer(null, variantAlleleTrigger);
  }

  public void testReference() {
    // homozygous identity
    final Variant v = createVariant("GTCCTAACT", "GTCCTAACT", null);
    assertEquals(v, getDecomposer().decompose(v).get(0));

    // heterozygous identity
    final Variant v2 = createVariant("GTCCTAACT", "GTCCTAACT", "GTCCTAACT");
    assertEquals(v2, getDecomposer().decompose(v2).get(0));
  }

  public void testHomozygousTrim() {
    // homozygous different lengths
    final Variant v3 = createVariant("GTCCTAACT", "GTCCTAAT", null);
    final List<Variant> list = getDecomposer().decompose(v3);
    assertEquals(1, list.size());
    final Variant first = list.get(0);
    assertEquals("", first.getSample(0).getName());
    assertEquals(30.5, first.getSample(0).getPosterior());
    assertEquals(false, first.getSample(0).isIdentity());
    final VariantLocus locus = first.getLocus();
    assertEquals("C", locus.getRefNts());
    assertEquals("chr", locus.getSequenceName());
    assertEquals(7, locus.getStart());
    assertEquals(8, locus.getEnd());
    assertEquals('A', locus.getPreviousRefNt());
  }

  public void testHeterozygousTrim() {
    // heterozygous different lengths
    final Variant v4 = createVariant("GTCCTAACT", "GTCCTAA", "GTCCTAACTAA");
    final List<Variant> list = getDecomposer().decompose(v4);
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
    //should not trim identity calls
    final Variant v = createVariant("ACC", "ACC", null);
    assertEquals(v, SimpleDecomposer.trim(v, TRIGGER));

    final Variant v2 = createVariant("ACC", "ACC", "ACC");
    assertEquals(v2, SimpleDecomposer.trim(v2, TRIGGER));
  }

  public void testTrimInsert() {
    //homozygous
    //                                               0123456789
    final Variant vHomozygousInsert = createVariant("ACGTCTGTCT", "ACGTTTCTGTCT", null);
    final Variant trim = SimpleDecomposer.trim(vHomozygousInsert, TRIGGER);
    assertEquals('G', trim.getLocus().getPreviousRefNt());
    assertEquals("", trim.getLocus().getRefNts());
    assertEquals(3, trim.getLocus().getStart());
    assertEquals(3, trim.getLocus().getEnd());
    assertEquals("TT", trim.getSample(0).getName());

    //heterozygous
    //                                               0123456789
    final Variant vHeterozyygousInsert = createVariant("ACGTCTGTCT", "ACGTTTCTGTCT", "ACGTCCCTGTCT");
    final Variant trim2 = SimpleDecomposer.trim(vHeterozyygousInsert, TRIGGER);
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
    final Variant trim = SimpleDecomposer.trim(vHomozygousIndel, TRIGGER);
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
    final Variant trim2 = SimpleDecomposer.trim(vHeterozyygousIndel, TRIGGER);
    assertEquals('C', trim2.getLocus().getPreviousRefNt());
    assertEquals("GTCTATC", trim2.getLocus().getRefNts());
    assertEquals(2, trim2.getLocus().getStart());
    assertEquals(9, trim2.getLocus().getEnd());
    assertEquals("CTATTTTC:GTCTAG", trim2.getSample(0).getName());
  }

  public void testVariantAlleleTrimOnly() {
    // There was a bug where the common descriptions weren't combined during trimming. This resulted in a VA = ref if
    // the third description outnumbered the variant here.
    final String name = "TAA:AAA";
    final boolean isIdentity = false;
    final String ref = "AAA";
    final String cat1 = "TAA";
    final String cat2 = "AAC";
    final VariantSample sample = VariantOutputVcfFormatterTest.createSample(Ploidy.DIPLOID, name, isIdentity, 30.5, VariantSample.DeNovoStatus.UNSPECIFIED, null);
    final String[] alleles = isIdentity ? new String[] {ref} : new String[] {ref, cat1, cat2};
    final DescriptionCommon description = new DescriptionCommon(alleles);
    final StatisticsComplex stats = new StatisticsComplex(description, 1);

    increment(description, stats, 0, 8);
    increment(description, stats, 1, 4);
    increment(description, stats, 2, 5);

    sample.setStats(stats);
    final Variant v = new Variant(new VariantLocus("chr", 0, ref.length(), ref, 'N'), sample);
    final VariantSample sample1 = v.getSample(0);
    sample1.setVariantAllele("TAA");

    final Variant split = SimpleDecomposer.trim(v, new VariantAlleleTrigger(1, 0.1));
    assertEquals("T", split.getSample(0).getVariantAllele());
  }
}
