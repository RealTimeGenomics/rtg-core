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

import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.rtg.reference.Ploidy;
import com.rtg.relation.ChildFamilyLookup;
import com.rtg.relation.Family;
import com.rtg.relation.LineageLookup;
import com.rtg.util.StringUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.variant.Variant;
import com.rtg.variant.VariantLocus;
import com.rtg.variant.VariantSample;
import com.rtg.variant.bayes.Description;
import com.rtg.variant.bayes.multisample.VariantAlleleTrigger;
import com.rtg.variant.bayes.snp.DescriptionCommon;
import com.rtg.variant.bayes.snp.EvidenceQ;
import com.rtg.variant.bayes.snp.StatisticsSnp;
import com.rtg.variant.format.VariantOutputVcfFormatterTest;
import com.rtg.variant.util.VariantUtils;

import junit.framework.TestCase;

/**
 */
public abstract class AbstractDecomposerTest extends TestCase {

  static Variant createVariant(final String ref, final String allele1, final String allele2) {
    final boolean hetero = allele2 != null;
    final boolean isIdentity = ref.equals(allele1) && (!hetero || ref.equals(allele2));
    final String name = allele1 + (hetero ? ":" + allele2 : "");
    final VariantSample sample = VariantOutputVcfFormatterTest.createSample(Ploidy.DIPLOID, name, isIdentity, 30.5, VariantSample.DeNovoStatus.UNSPECIFIED, null);
    final String[] alleles = isIdentity ? new String[] {ref} : hetero ? new String[] {ref, allele1, allele2} : new String[] {ref, allele1};
    sample.setStats(new StatisticsComplex(new DescriptionCommon(alleles), 1));
    return new Variant(new VariantLocus("chr", 0, ref.length(), ref, 'N'), sample);
  }

  static VariantSample getVariantSample(Ploidy diploid, String sequence, boolean b, Double score, VariantSample.DeNovoStatus isDeNovo, Double deNovoPosterior, Description desc) {
    final VariantSample vs = VariantOutputVcfFormatterTest.createSample(diploid, sequence, b, score, isDeNovo, deNovoPosterior);
    vs.setStats(new StatisticsSnp(desc));
    return vs;
  }

  static DescriptionCommon getDescription(final String... names) {
    final HashSet<String> allelesSet = new HashSet<>();
    for (final String name : names) {
      allelesSet.addAll(Arrays.asList(StringUtils.split(name, VariantUtils.COLON)));
    }
    return new DescriptionCommon(allelesSet.toArray(new String[allelesSet.size()]));
  }

  static void checkDenovoCorrect(List<Variant> variants, int[] expectedPositions, String[] expectedParentA, String[] expectedParentB, String[] expectedChild, VariantSample.DeNovoStatus[] expectedDenovo) {
    for (int i = 0; i < variants.size(); ++i) {
      final Variant res = variants.get(i);
      assertEquals(expectedPositions[i], res.getLocus().getStart());
      assertEquals(res.getLocus().getStart() + 1, res.getLocus().getEnd());
      if (expectedParentA[i] != null) {
        assertEquals(expectedParentA[i], res.getSample(0).getName());
      } else {
        assertNull(res.getSample(0));
      }
      if (expectedParentB[i] != null) {
        assertEquals(expectedParentB[i], res.getSample(1).getName());
      } else {
        assertNull(res.getSample(1));
      }
      assertEquals(expectedChild[i], res.getSample(2).getName());
      assertEquals(expectedDenovo[i], res.getSample(2).isDeNovo());
    }
  }

  static void increment(DescriptionCommon description, StatisticsComplex stats, int read, int count) {
    for (int i = 0; i < count; ++i) {
      stats.increment(new EvidenceQ(description, read, 2, 2, 0.1, 0.1, true, true, true, false), 0);
    }
  }

  public abstract Decomposer getDecomposer();

  public abstract Decomposer getDecomposer(DenovoChecker denovoChecker);

  public abstract Decomposer getDecomposer(VariantAlleleTrigger variantAlleleTrigger);

  public void testSplitHomozygous() {
    final String a = "GTCCTAACT";
    final String b = "GCTCTCACG";
    final Variant v = createVariant(a, b, null);

    final Map<Set<String>, Double> likelihoods = new HashMap<>();
    likelihoods.put(Collections.singleton(a), Math.log(0.1));
    likelihoods.put(VariantSample.pairSet(a, b), Math.log(0.5));
    likelihoods.put(Collections.singleton(b), Math.log(0.4));
    v.getSample(0).setGenotypeLikelihoods(likelihoods);

    final List<Variant> list = getDecomposer().decompose(v);
    assertEquals(3, list.size());
    final Variant first = list.get(0);
    assertEquals("CT", first.getSample(0).getName());
    assertEquals(30.5, first.getSample(0).getPosterior());
    assertEquals(false, first.getSample(0).isIdentity());
    assertEquals("TC", first.getLocus().getRefNts());
    assertEquals("chr", first.getLocus().getSequenceName());
    assertEquals(1, first.getLocus().getStart());
    assertEquals(3, first.getLocus().getEnd());
    assertEquals('G', first.getLocus().getPreviousRefNt());
    assertEquals(Math.log(0.4), first.getSample(0).getGenotypeLikelihoods().get(Collections.singleton("CT")));

    final Variant second = list.get(1);
    assertEquals("C", second.getSample(0).getName());
    assertEquals(30.5, second.getSample(0).getPosterior());
    assertEquals(false, second.getSample(0).isIdentity());
    assertEquals("A", second.getLocus().getRefNts());
    assertEquals("chr", second.getLocus().getSequenceName());
    assertEquals(5, second.getLocus().getStart());
    assertEquals(6, second.getLocus().getEnd());
    assertEquals('T', second.getLocus().getPreviousRefNt());
    assertEquals(Math.log(0.5), second.getSample(0).getGenotypeLikelihoods().get(VariantSample.pairSet("C", "A")));

    final Variant third = list.get(2);
    assertEquals("G", third.getSample(0).getName());
    assertEquals(30.5, third.getSample(0).getPosterior());
    assertEquals(false, third.getSample(0).isIdentity());
    assertEquals("T", third.getLocus().getRefNts());
    assertEquals("chr", third.getLocus().getSequenceName());
    assertEquals(8, third.getLocus().getStart());
    assertEquals(9, third.getLocus().getEnd());
    assertEquals('C', third.getLocus().getPreviousRefNt());
    assertEquals(Math.log(0.1), third.getSample(0).getGenotypeLikelihoods().get(Collections.singleton("T")));
  }

  public void testSplitHeterozygous() {
    final String a = "GTCCTAACT";
    final String b = "GCTCTCACG";
    final String c = "GTCCTAACG";
    final Variant v = createVariant(a, b, c);

    final Map<Set<String>, Double> likelihoods = new HashMap<>();
    likelihoods.put(Collections.singleton(a), Math.log(0.1));
    likelihoods.put(Collections.singleton(b), Math.log(0.12));
    likelihoods.put(Collections.singleton(c), Math.log(0.105));
    likelihoods.put(VariantSample.pairSet(a, b), Math.log(0.15));
    likelihoods.put(VariantSample.pairSet(a, c), Math.log(0.18));
    likelihoods.put(VariantSample.pairSet(b, c), Math.log(0.345));
    v.getSample(0).setGenotypeLikelihoods(likelihoods);

    final List<Variant> list = getDecomposer().decompose(v);
    assertEquals(3, list.size());
    final Variant first = list.get(0);
    assertEquals("CT:TC", first.getSample(0).getName());
    assertEquals(30.5, first.getSample(0).getPosterior());
    assertEquals(false, first.getSample(0).isIdentity());
    assertEquals("TC", first.getLocus().getRefNts());
    assertEquals("chr", first.getLocus().getSequenceName());
    assertEquals(1, first.getLocus().getStart());
    assertEquals(3, first.getLocus().getEnd());
    assertEquals('G', first.getLocus().getPreviousRefNt());
    assertEquals(0, first.getSplitId());
    assertEquals(Math.log(0.385), first.getSample(0).getGenotypeLikelihoods().get(Collections.singleton("TC")), 0.0001);
    assertEquals(Math.log(0.12), first.getSample(0).getGenotypeLikelihoods().get(Collections.singleton("CT")), 0.0001);
    assertEquals(Math.log(0.495), first.getSample(0).getGenotypeLikelihoods().get(VariantSample.pairSet("CT", "TC")), 0.0001);

    final Variant second = list.get(1);
    assertEquals("C:A", second.getSample(0).getName());
    assertEquals(30.5, second.getSample(0).getPosterior());
    assertEquals(false, second.getSample(0).isIdentity());
    assertEquals("A", second.getLocus().getRefNts());
    assertEquals("chr", second.getLocus().getSequenceName());
    assertEquals(5, second.getLocus().getStart());
    assertEquals(6, second.getLocus().getEnd());
    assertEquals('T', second.getLocus().getPreviousRefNt());
    assertEquals(1, second.getSplitId());

    final Variant third = list.get(2);
    assertEquals("G:G", third.getSample(0).getName());
    assertEquals(30.5, third.getSample(0).getPosterior());
    assertEquals(false, third.getSample(0).isIdentity());
    assertEquals("T", third.getLocus().getRefNts());
    assertEquals("chr", third.getLocus().getSequenceName());
    assertEquals(8, third.getLocus().getStart());
    assertEquals(9, third.getLocus().getEnd());
    assertEquals('C', third.getLocus().getPreviousRefNt());
    assertEquals(2, third.getSplitId());
  }

  public void testSplitTrim() {
    final VariantLocus locus = new VariantLocus("blah", 1, 8, "ACTACAG", '?');
    final VariantSample vs = getVariantSample(Ploidy.DIPLOID, "AGTACAG:AGTACAG", false, null, VariantSample.DeNovoStatus.NOT_DE_NOVO, 0.0, getDescription(locus.getRefNts(), "AGTACAG:AGTACAG"));
    vs.setCoverage(10);
    vs.setCoverageCorrection(1.0);
    final Variant v = new Variant(locus, vs);
    v.setNonIdentityPosterior(2.3);
    final List<Variant> list = getDecomposer().decompose(v);
    assertEquals(1, list.size());
    final Variant variant = list.get(0);
    assertEquals('A', variant.getLocus().getPreviousRefNt());
    assertEquals(2, variant.getLocus().getStart());
    assertEquals(3, variant.getLocus().getEnd());
  }

  public void testSplitWithSomaticCause() {
    final VariantLocus locus = new VariantLocus("blah", 1, 4, "AGT", '?');
    final VariantSample normal = getVariantSample(Ploidy.DIPLOID, "CGA:CGA", false, null, VariantSample.DeNovoStatus.UNSPECIFIED, null, getDescription(locus.getRefNts(), "CGA:CGT"));
    normal.setCoverage(10);
    normal.setCoverageCorrection(1.0);
    final VariantSample cancer = getVariantSample(Ploidy.DIPLOID, "CGA:CGT", false, null, VariantSample.DeNovoStatus.IS_DE_NOVO, 10.0, getDescription(locus.getRefNts(), "CGA:CGT"));
    cancer.setCoverage(10);
    cancer.setCoverageCorrection(1.0);
    final Variant v = new Variant(locus, normal, cancer);
    v.setNonIdentityPosterior(2.3);
    final List<Variant> list = getDecomposer(new LineageDenovoChecker(new LineageLookup(-1, 0))).decompose(v);
    assertEquals(2, list.size());
    assertEquals(VariantSample.DeNovoStatus.NOT_DE_NOVO, list.get(0).getSample(1).isDeNovo());
    assertEquals(10.0, list.get(1).getSample(1).getDeNovoPosterior());
    assertEquals(VariantSample.DeNovoStatus.IS_DE_NOVO, list.get(1).getSample(1).isDeNovo());
  }

  public void testOverCoverageGetsSet() {
    Diagnostic.setLogStream();
    final VariantLocus locus = new VariantLocus("blah", 0, 5, "cgtgt", 'N');
    final VariantSample vs = getVariantSample(Ploidy.DIPLOID, "cgt", false, null, VariantSample.DeNovoStatus.NOT_DE_NOVO, 0.0, getDescription(locus.getRefNts(), "cgt"));
    vs.setCoverage(3);
    vs.setCoverageCorrection(0.000);
    final Variant v = new Variant(locus, vs);
    v.addFilter(Variant.VariantFilter.COVERAGE);
    final List<Variant> list = getDecomposer().decompose(v);
    assertEquals(1, list.size());
    final Variant variant = list.get(0);
    assertTrue(variant.isFiltered(Variant.VariantFilter.COVERAGE));
  }

  public void testPrevNt() {
    final VariantLocus locus = new VariantLocus("blah", 0, 3, "CGT", 'N');
    final VariantSample vs = getVariantSample(Ploidy.DIPLOID, "CGTGT", false, null, VariantSample.DeNovoStatus.NOT_DE_NOVO, 0.0, getDescription("CGTGT", locus.getRefNts()));
    final Variant v = new Variant(locus, vs);
    final List<Variant> list = getDecomposer().decompose(v);
    assertEquals(1, list.size());
    final VariantLocus newLocus = list.get(0).getLocus();
    assertEquals("GT", list.get(0).getSample(0).getName());
    assertEquals('C', newLocus.getPreviousRefNt());
    assertEquals("", newLocus.getRefNts());
    assertEquals(1, newLocus.getStart());
  }

  public void testTrimSimpleSnps() {
    // homozygous
    final Variant trim = getDecomposer().decompose(createVariant("ACGTCTGTCT", "ACCTCTGTCT", null)).get(0);
    assertEquals('C', trim.getLocus().getPreviousRefNt());
    assertEquals("G", trim.getLocus().getRefNts());
    assertEquals(2, trim.getLocus().getStart());
    assertEquals(3, trim.getLocus().getEnd());
    assertEquals("C", trim.getSample(0).getName());

    // heterozygous
    final Variant trim2 = getDecomposer().decompose(createVariant("ACGTCTGTCT", "ACCTCTGTCT", "ACGTCTGTCT")).get(0);
    assertEquals('C', trim2.getLocus().getPreviousRefNt());
    assertEquals("G", trim2.getLocus().getRefNts());
    assertEquals(2, trim2.getLocus().getStart());
    assertEquals(3, trim2.getLocus().getEnd());
    assertEquals("C:G", trim2.getSample(0).getName());
  }

  public void testTrimDelete() {
    // homozygous delete
    final Variant trim = getDecomposer().decompose(createVariant("ACGTCTGTCT", "ACGTCTGT", null)).get(0);
    assertEquals('G', trim.getLocus().getPreviousRefNt());
    assertEquals("TC", trim.getLocus().getRefNts());
    assertEquals(7, trim.getLocus().getStart());
    assertEquals(9, trim.getLocus().getEnd());
    assertEquals("", trim.getSample(0).getName());

    // heterozygous deletes
    final Variant trim2 = getDecomposer().decompose(createVariant("ACGTCTGTCT", "ACGTCTGT", "ACGTCTG")).get(0);
    assertEquals('G', trim2.getLocus().getPreviousRefNt());
    assertEquals("TCT", trim2.getLocus().getRefNts());
    assertEquals(7, trim2.getLocus().getStart());
    assertEquals(10, trim2.getLocus().getEnd());
    assertEquals("T:", trim2.getSample(0).getName());

    // heterozygous deletes
    final Variant trim3 = getDecomposer().decompose(createVariant("ACGTCTGTCT", "ACGTCTG", "ACGTCTGT")).get(0);
    assertEquals('G', trim3.getLocus().getPreviousRefNt());
    assertEquals("TCT", trim3.getLocus().getRefNts());
    assertEquals(":T", trim3.getSample(0).getName());
  }

  public void testTrimInsertHomo() {
    final Variant trim = getDecomposer().decompose(createVariant("ACGTCTGTCT", "ACGTTTCTGTCT", null)).get(0);
    assertEquals("TT", trim.getSample(0).getName());
    assertEquals('G', trim.getLocus().getPreviousRefNt());
    assertEquals("", trim.getLocus().getRefNts());
    assertEquals(3, trim.getLocus().getStart());
    assertEquals(3, trim.getLocus().getEnd());
  }

  public void testTrimInsertHetero() {
    final List<Variant> decomposed = getDecomposer().decompose(createVariant("ACGTCTGTCT", "ACGTTTCTGTCT", "ACGTCCCTGTCT"));
    final Variant trim2 = decomposed.get(0);
    assertEquals("TT:CC", trim2.getSample(0).getName());
    assertEquals('T', trim2.getLocus().getPreviousRefNt());
    assertEquals("", trim2.getLocus().getRefNts());
    assertEquals(4, trim2.getLocus().getStart());
    assertEquals(4, trim2.getLocus().getEnd());
  }

  public void testTrimComplexMultisample() {
    final VariantLocus locus = new VariantLocus("chr", 0, "ACGTCCCCT".length(), "ACGTCCCCT", 'N');
    final DescriptionCommon description = getDescription(locus.getRefNts(), "ACGTCTGTCT", "ACGTCTGTCT:ACGTTTT", "ACGTCTCT");
    final VariantSample sample1 = getVariantSample(Ploidy.DIPLOID, "ACGTCTGTCT", false, 30.5, VariantSample.DeNovoStatus.UNSPECIFIED, null, description);
    final VariantSample sample2 = getVariantSample(Ploidy.DIPLOID, "ACGTCTGTCT:ACGTTTT", false, 30.5, VariantSample.DeNovoStatus.UNSPECIFIED, null, description);
    final VariantSample sample3 = getVariantSample(Ploidy.DIPLOID, "ACGTCTCT", false, 30.5, VariantSample.DeNovoStatus.UNSPECIFIED, null, description);
    final Variant v = new Variant(locus, sample1, sample2, sample3);
    final Variant trim = getDecomposer().decompose(v).get(0);
    assertEquals('T', trim.getLocus().getPreviousRefNt());
    assertEquals("CCCC", trim.getLocus().getRefNts());
    assertEquals(4, trim.getLocus().getStart());
    assertEquals(8, trim.getLocus().getEnd());
    assertEquals("CTGTC", trim.getSample(0).getName());
    assertEquals("CTGTC:TT", trim.getSample(1).getName());
    assertEquals("CTC", trim.getSample(2).getName());
  }

  public void testTrimComplexMultisampleWithNulls() {
    final VariantLocus locus = new VariantLocus("chr", 0, "ACGTCCCCT".length(), "ACGTCCCCT", 'N');
    final VariantSample sample1 = getVariantSample(Ploidy.DIPLOID, "ACGTCTGTCT", false, 30.5, VariantSample.DeNovoStatus.UNSPECIFIED, null, getDescription(locus.getRefNts(), "ACGTCTGTCT"));
    final Variant v = new Variant(locus, sample1, null, null);
    final Variant trim = getDecomposer().decompose(v).get(0);
    assertEquals('C', trim.getLocus().getPreviousRefNt());
    assertEquals("CC", trim.getLocus().getRefNts());
    assertEquals(5, trim.getLocus().getStart());
    assertEquals(7, trim.getLocus().getEnd());
    assertEquals("TGT", trim.getSample(0).getName());
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

  public void testDenovoCorrect() {
    final VariantLocus locus = new VariantLocus("chr", 0, "TTTT".length(), "TTTT", 'N');
    final DescriptionCommon description = AligningDecomposerTest.getDescription(locus.getRefNts(), "TTTT:TTTC", "TTTT:TCTC");
    final VariantSample vsDad = getVariantSample(Ploidy.DIPLOID, "TTTT:TTTC", false, 5.0, VariantSample.DeNovoStatus.NOT_DE_NOVO, 0.0, description);
    final VariantSample vsMum = getVariantSample(Ploidy.DIPLOID, "TTTT:TTTC", false, 5.0, VariantSample.DeNovoStatus.NOT_DE_NOVO, 0.0, description);
    final VariantSample vsJimmy = getVariantSample(Ploidy.DIPLOID, "TTTT:TCTC", false, 5.0, VariantSample.DeNovoStatus.IS_DE_NOVO, null, description);
    final Variant v = new Variant(locus, vsDad, vsMum, vsJimmy);
    final List<Variant> variants = getDecomposer(new MendelianDenovoChecker(new ChildFamilyLookup(3, new Family("father", "mother", "child")))).decompose(v);
    final int[] expectedPositions = {1, 3};
    final String[] expectedParentA = {"T:T", "T:C"};
    final String[] expectedParentB = {"T:T", "T:C"};
    final String[] expectedChild = {"T:C", "T:C"};
    final VariantSample.DeNovoStatus[] expectedDenovo = new VariantSample.DeNovoStatus[] {VariantSample.DeNovoStatus.IS_DE_NOVO, VariantSample.DeNovoStatus.NOT_DE_NOVO};
    checkDenovoCorrect(variants, expectedPositions, expectedParentA, expectedParentB, expectedChild, expectedDenovo);
  }

  public void testDenovoCorrectHaploid() {
    final VariantLocus locus = new VariantLocus("chr", 0, "TTTT".length(), "TTTT", 'N');
    final DescriptionCommon description = getDescription(locus.getRefNts(), "TTTT", "TTTT:TTTC", "TCTC");
    final VariantSample vsDad = getVariantSample(Ploidy.HAPLOID, "TTTT", false, 5.0, VariantSample.DeNovoStatus.NOT_DE_NOVO, 0.0, description);
    final VariantSample vsMum = getVariantSample(Ploidy.DIPLOID, "TTTT:TTTC", false, 5.0, VariantSample.DeNovoStatus.NOT_DE_NOVO, 0.0, description);
    final VariantSample vsJimmy = getVariantSample(Ploidy.HAPLOID, "TCTC", false, 5.0, VariantSample.DeNovoStatus.IS_DE_NOVO, null, description);
    final Variant v = new Variant(locus, vsDad, vsMum, vsJimmy);
    final List<Variant> variants = getDecomposer(new MendelianDenovoChecker(new ChildFamilyLookup(3, new Family("father", "mother", "child")))).decompose(v);
    final int[] expectedPositions = {1, 3};
    final String[] expectedParentA = {"T", "T"};
    final String[] expectedParentB = {"T:T", "T:C"};
    final String[] expectedChild = {"C", "C"};
    final VariantSample.DeNovoStatus[] expectedDenovo = new VariantSample.DeNovoStatus[] {VariantSample.DeNovoStatus.IS_DE_NOVO, VariantSample.DeNovoStatus.NOT_DE_NOVO};
    checkDenovoCorrect(variants, expectedPositions, expectedParentA, expectedParentB, expectedChild, expectedDenovo);
  }

  public void testDenovoCorrectHaploid2() {
    final VariantLocus locus = new VariantLocus("chr", 0, "TTTT".length(), "TTTT", 'N');
    final DescriptionCommon description = getDescription(locus.getRefNts(), "TTTC", "TCTC");
    final VariantSample vsDad = getVariantSample(Ploidy.HAPLOID, "TTTC", false, 5.0, VariantSample.DeNovoStatus.NOT_DE_NOVO, 0.0, description);
    final VariantSample vsJimmy = getVariantSample(Ploidy.HAPLOID, "TCTC", false, 5.0, VariantSample.DeNovoStatus.IS_DE_NOVO, null, description);
    final Variant v = new Variant(locus, vsDad, null, vsJimmy);
    final List<Variant> variants = getDecomposer(new MendelianDenovoChecker(new ChildFamilyLookup(3, new Family("father", "mother", "child")))).decompose(v);
    final int[] expectedPositions = {1, 3};
    final String[] expectedParentA = {"T", "C"};
    final String[] expectedParentB = {null, null};
    final String[] expectedChild = {"C", "C"};
    final VariantSample.DeNovoStatus[] expectedDenovo = new VariantSample.DeNovoStatus[] {VariantSample.DeNovoStatus.IS_DE_NOVO, VariantSample.DeNovoStatus.NOT_DE_NOVO};
    checkDenovoCorrect(variants, expectedPositions, expectedParentA, expectedParentB, expectedChild, expectedDenovo);
  }

  public void testVariantAlleleNearHomozygous() {
    final String name = "TAG:TAA";
    final String ref = "AAA";
    final String cat1 = "TAA";
    final String cat2 = "TAG";
    final VariantSample sample = VariantOutputVcfFormatterTest.createSample(Ploidy.DIPLOID, name, false, 30.5, VariantSample.DeNovoStatus.UNSPECIFIED, null);
    final DescriptionCommon description = new DescriptionCommon(ref, cat1, cat2);
    final StatisticsComplex stats = new StatisticsComplex(description, 1);

    increment(description, stats, 0, 4);
    increment(description, stats, 1, 3);
    increment(description, stats, 2, 2);

    sample.setStats(stats);
    final Variant v = new Variant(new VariantLocus("chr", 0, ref.length(), ref, 'N'), sample);
    final VariantSample sample1 = v.getSample(0);
    sample1.setVariantAllele("TAG");

    final List<Variant> splits = getDecomposer(new VariantAlleleTrigger(1, 0.1)).decompose(v);
    final VariantSample s0 = splits.get(0).getSample(0);
    assertEquals("T:T", s0.getName());
    assertEquals("T", s0.getVariantAllele());
    final VariantSample s1 = splits.get(1).getSample(0);
    assertEquals("G:A", s1.getName());
    assertEquals("G", s1.getVariantAllele());
  }

  public void testVariantAlleleCrossedSteams() {
    final String name = "AAG:TAA";
    final boolean isIdentity = false;
    final String ref = "AAA";
    final String cat1 = "AAG";
    final String cat2 = "TAA";
    final VariantSample sample = VariantOutputVcfFormatterTest.createSample(Ploidy.DIPLOID, name, isIdentity, 30.5, VariantSample.DeNovoStatus.UNSPECIFIED, null);
    final String[] alleles = isIdentity ? new String[] {ref} : new String[] {ref, cat1, cat2};
    final DescriptionCommon description = new DescriptionCommon(alleles);
    final StatisticsComplex stats = new StatisticsComplex(description, 1);

    increment(description, stats, 0, 4);
    increment(description, stats, 1, 3);
    increment(description, stats, 2, 2);

    sample.setStats(stats);
    final Variant v = new Variant(new VariantLocus("chr", 0, ref.length(), ref, 'N'), sample);
    final VariantSample sample1 = v.getSample(0);
    sample1.setVariantAllele("TAA");

    final List<Variant> splits = getDecomposer(new VariantAlleleTrigger(1, 0.1)).decompose(v);
    assertEquals("T", splits.get(0).getSample(0).getVariantAllele());
    assertEquals("G", splits.get(1).getSample(0).getVariantAllele());
  }
}
