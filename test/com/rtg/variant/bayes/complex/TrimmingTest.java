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
import com.rtg.util.StringUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.variant.StaticThreshold;
import com.rtg.variant.Variant;
import com.rtg.variant.Variant.VariantFilter;
import com.rtg.variant.VariantLocus;
import com.rtg.variant.VariantParams;
import com.rtg.variant.VariantSample;
import com.rtg.variant.bayes.Description;
import com.rtg.variant.bayes.snp.DescriptionCommon;
import com.rtg.variant.bayes.snp.StatisticsSnp;
import com.rtg.variant.format.VariantOutputVcfFormatterTest;

import junit.framework.TestCase;

/**
 */
public class TrimmingTest extends TestCase {

  public void testSplitHomozygous() {
    final String a = "GTCCTAACT";
    final String b = "GCTCTCACG";
    final Variant v = createVariant(a, b, null);

    final Map<Set<String>, Double> likelihoods = new HashMap<>();
    likelihoods.put(Collections.singleton(a), Math.log(0.1));
    likelihoods.put(VariantSample.pairSet(a, b), Math.log(0.5));
    likelihoods.put(Collections.singleton(b), Math.log(0.4));
    v.getSample(0).setGenotypeLikelihoods(likelihoods);

    final List<Variant> list = Trimming.split(v, null);
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

    final List<Variant> list = Trimming.split(v, null);
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

  public void testNonSplit() {
    //homozygous identity
    final Variant v = createVariant("GTCCTAACT", "GTCCTAACT", null);
    final List<Variant> list = Trimming.split(v, null);
    assertEquals(v, list.get(0));

    //heterozygous identity
    final Variant v2 = createVariant("GTCCTAACT", "GTCCTAACT", "GTCCTAACT");
    final List<Variant> list2 = Trimming.split(v2, null);
    assertEquals(v2, list2.get(0));

    //homozygous different lengths
    final Variant v3 = createVariant("GTCCTAACT", "GTCCTAAT", null);
    final List<Variant> list3 = Trimming.split(v3, null);
    assertEquals(v3, list3.get(0));

    //heterozygous different lengths
    final Variant v4 = createVariant("GTCCTAACT", "GTCCTAA", "GTCCTAACTAA");
    final List<Variant> list4 = Trimming.split(v4, null);
    assertEquals(v4, list4.get(0));
  }

  public void testPrevNt() {
    Diagnostic.setLogStream();
    final VariantParams p = VariantParams.builder().create();
    final VariantLocus locus = new VariantLocus("blah", 0, 3, "cgt", 'N');
    final VariantSample vs = getVariantSample(Ploidy.DIPLOID, "cgtgt", false, null, VariantSample.DeNovoStatus.NOT_DE_NOVO, 0.0, getDescription("cgtgt", locus.getRefNts()));
    final Variant v = new Variant(locus, vs);
    final List<Variant> list = Trimming.trimSplit(p, v, null);
    assertEquals(1, list.size());
    final VariantLocus newLocus = list.get(0).getLocus();
    assertEquals("gt", list.get(0).getSample(0).getName());
    assertEquals('c', newLocus.getPreviousRefNt());
    assertEquals("", newLocus.getRefNts());
    assertEquals(1, newLocus.getStart());
  }

  static VariantSample getVariantSample(Ploidy diploid, String cgtgt, boolean b, Double score, VariantSample.DeNovoStatus isDeNovo, Double deNovoPosterior, Description desc) {
    final VariantSample vs = VariantOutputVcfFormatterTest.createSample(diploid, cgtgt, b, score, isDeNovo, deNovoPosterior);
    vs.setStats(new StatisticsSnp(desc));
    return vs;
  }

  static DescriptionCommon getDescription(String... names) {
    final HashSet<String> allelesSet = new HashSet<>();
    for (String name : names) {
      allelesSet.addAll(Arrays.asList(StringUtils.split(name, ':')));
    }
    return new DescriptionCommon(allelesSet.toArray(new String[allelesSet.size()]));
  }

  public void testSplitTrim() {
    final VariantParams p = VariantParams.builder().outputNonSnps(false).create();
    final VariantLocus locus = new VariantLocus("blah", 1, 8, "ACTACAG", '?');
    final VariantSample vs = getVariantSample(Ploidy.DIPLOID, "AGTACAG:AGTACAG", false, null, VariantSample.DeNovoStatus.NOT_DE_NOVO, 0.0, getDescription(locus.getRefNts(), "AGTACAG:AGTACAG"));
    vs.setCoverage(10);
    vs.setCoverageCorrection(1.0);
    final Variant v = new Variant(locus, vs);
    v.setNonIdentityPosterior(2.3);
    final List<Variant> list = Trimming.trimSplit(p, v, null);
    assertEquals(1, list.size());
    final Variant variant = list.get(0);
    assertEquals('A', variant.getLocus().getPreviousRefNt());
    assertEquals(2, variant.getLocus().getStart());
    assertEquals(3, variant.getLocus().getEnd());
  }

  public void testSplitWithSomaticCause() {
    final VariantParams p = VariantParams.builder().outputNonSnps(false).create();
    final VariantLocus locus = new VariantLocus("blah", 1, 4, "AGT", '?');
    final VariantSample vs = getVariantSample(Ploidy.DIPLOID, "CGA:CGT", false, null, VariantSample.DeNovoStatus.NOT_DE_NOVO, 0.0, getDescription(locus.getRefNts(), "CGA:CGT"));
    vs.setCoverage(10);
    vs.setCoverageCorrection(1.0);
    final Variant v = new Variant(locus, vs);
    v.setNonIdentityPosterior(2.3);
    v.setPossibleCause("CGA:CGT");
    final List<Variant> list = Trimming.trimSplit(p, v, null);
    assertEquals(2, list.size());
    assertNull(list.get(0).getPossibleCause());
    assertEquals("A:T", list.get(1).getPossibleCause());
  }

  public void testOverCoverageGetsSet() {
    Diagnostic.setLogStream();
    final VariantParams p = VariantParams.builder().maxCoverageFilter(new StaticThreshold(2)).create();
    final VariantLocus locus = new VariantLocus("blah", 0, 5, "cgtgt", 'N');
    final VariantSample vs = getVariantSample(Ploidy.DIPLOID, "cgt", false, null, VariantSample.DeNovoStatus.NOT_DE_NOVO, 0.0, getDescription(locus.getRefNts(), "cgt"));
    vs.setCoverage(3);
    vs.setCoverageCorrection(0.000);
    final Variant v = new Variant(locus, vs);
    v.addFilter(VariantFilter.COVERAGE);
    final List<Variant> list = Trimming.trimSplit(p, v, null);
    assertEquals(1, list.size());
    final Variant varianceCall = list.get(0);
    assertTrue(varianceCall.isFiltered(VariantFilter.COVERAGE));
  }

  private static Variant createVariant(String ref, String cat1, String cat2) {
    final boolean hetero = cat2 != null;
    final boolean isIdentity = ref.equals(cat1) && (!hetero || ref.equals(cat2));
    final String name = cat1 + (hetero ? ":" + cat2 : "");
    final VariantSample sample = VariantOutputVcfFormatterTest.createSample(Ploidy.DIPLOID, name, isIdentity, 30.5, VariantSample.DeNovoStatus.UNSPECIFIED, null);
    final String[] alleles = isIdentity ? new String[] {ref} : hetero ? new String[] {ref, cat1, cat2} : new String[] {ref, cat1};
    sample.setStats(new StatisticsComplex(new DescriptionCommon(alleles), 1));
    return new Variant(new VariantLocus("chr", 0, ref.length(), ref, 'N'), sample);
  }

  public void testTrimIdentity() {
    //should not trim identity calls
    final Variant v = createVariant("ACC", "ACC", null);
    assertEquals(v, Trimming.trim(v));

    final Variant v2 = createVariant("ACC", "ACC", "ACC");
    assertEquals(v2, Trimming.trim(v2));
  }

  public void testTrimSimpleSnps() {
    //homozygous
    //                               0123456789
    final Variant v = createVariant("ACGTCTGTCT", "ACCTCTGTCT", null);
    final Variant trim = Trimming.trim(v);
    assertEquals('C', trim.getLocus().getPreviousRefNt());
    assertEquals("G", trim.getLocus().getRefNts());
    assertEquals(2, trim.getLocus().getStart());
    assertEquals(3, trim.getLocus().getEnd());
    assertEquals("C", trim.getSample(0).getName());

    //heterozygous
    //                                0123456789
    final Variant v2 = createVariant("ACGTCTGTCT", "ACCTCTGTCT", "ACGTCTGTCT");
    final Variant trim2 = Trimming.trim(v2);
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
    final Variant trim = Trimming.trim(vHomozygousDelete);
    assertEquals('G', trim.getLocus().getPreviousRefNt());
    assertEquals("TC", trim.getLocus().getRefNts());
    assertEquals(7, trim.getLocus().getStart());
    assertEquals(9, trim.getLocus().getEnd());
    assertEquals("", trim.getSample(0).getName());

    //heterozygous deletes
    final Variant vHeterozygousDelete = createVariant("ACGTCTGTCT", "ACGTCTGT", "ACGTCTG");
    final Variant trim2 = Trimming.trim(vHeterozygousDelete);
    assertEquals('G', trim2.getLocus().getPreviousRefNt());
    assertEquals("TCT", trim2.getLocus().getRefNts());
    assertEquals(7, trim2.getLocus().getStart());
    assertEquals(10, trim2.getLocus().getEnd());
    assertEquals("T:", trim2.getSample(0).getName());

  //heterozygous deletes
    final Variant vHeterozygousDelete2 = createVariant("ACGTCTGTCT", "ACGTCTG", "ACGTCTGT");
    final Variant trim3 = Trimming.trim(vHeterozygousDelete2);
    assertEquals('G', trim3.getLocus().getPreviousRefNt());
    assertEquals("TCT", trim3.getLocus().getRefNts());
    assertEquals(":T", trim3.getSample(0).getName());
  }

  public void testTrimInsert() {
    //homozygous
    //                                               0123456789
    final Variant vHomozygousInsert = createVariant("ACGTCTGTCT", "ACGTTTCTGTCT", null);
    final Variant trim = Trimming.trim(vHomozygousInsert);
    assertEquals('G', trim.getLocus().getPreviousRefNt());
    assertEquals("", trim.getLocus().getRefNts());
    assertEquals(3, trim.getLocus().getStart());
    assertEquals(3, trim.getLocus().getEnd());
    assertEquals("TT", trim.getSample(0).getName());

    //heterozygous
    //                                               0123456789
    final Variant vHeterozyygousInsert = createVariant("ACGTCTGTCT", "ACGTTTCTGTCT", "ACGTCCCTGTCT");
    final Variant trim2 = Trimming.trim(vHeterozyygousInsert);
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
    final Variant trim = Trimming.trim(vHomozygousIndel);
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
    final Variant trim2 = Trimming.trim(vHeterozyygousIndel);
    assertEquals('C', trim2.getLocus().getPreviousRefNt());
    assertEquals("GTCTATC", trim2.getLocus().getRefNts());
    assertEquals(2, trim2.getLocus().getStart());
    assertEquals(9, trim2.getLocus().getEnd());
    assertEquals("CTATTTTC:GTCTAG", trim2.getSample(0).getName());
  }

  public void testTrimComplexMultisample() {
    final VariantLocus locus = new VariantLocus("chr", 0, "ACGTCCCCT".length(), "ACGTCCCCT", 'N');
    final DescriptionCommon description = getDescription(locus.getRefNts(), "ACGTCTGTCT", "ACGTCTGTCT:ACGTTTT", "ACGTCTCT");
    final VariantSample sample1 = getVariantSample(Ploidy.DIPLOID, "ACGTCTGTCT", false, 30.5, VariantSample.DeNovoStatus.UNSPECIFIED, null, description);
    final VariantSample sample2 = getVariantSample(Ploidy.DIPLOID, "ACGTCTGTCT:ACGTTTT", false, 30.5, VariantSample.DeNovoStatus.UNSPECIFIED, null, description);
    final VariantSample sample3 = getVariantSample(Ploidy.DIPLOID, "ACGTCTCT", false, 30.5, VariantSample.DeNovoStatus.UNSPECIFIED, null, description);
    final Variant v = new Variant(locus, sample1, sample2, sample3);
    final Variant trim = Trimming.trim(v);
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
    final Variant trim = Trimming.trim(v);
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
    final Variant trim = Trimming.trim(vHomozygousIndel);
    assertEquals('T', trim.getLocus().getPreviousRefNt());
    assertEquals("", trim.getLocus().getRefNts());
    assertEquals(1, trim.getLocus().getStart());
    assertEquals(1, trim.getLocus().getEnd());
    assertEquals("A", trim.getSample(0).getName());

    //heterozygous
    // AC|GTCTA----G|T
    final Variant vHeterozyygousIndel = createVariant("TA", "TAA", "TA");
    final Variant trim2 = Trimming.trim(vHeterozyygousIndel);
    assertEquals('T', trim2.getLocus().getPreviousRefNt());
    assertEquals("", trim2.getLocus().getRefNts());
    assertEquals(1, trim2.getLocus().getStart());
    assertEquals(1, trim2.getLocus().getEnd());
    assertEquals("A:", trim2.getSample(0).getName());
  }

  public void testDenovoCorrect() {
    final VariantLocus locus = new VariantLocus("chr", 0, "TTTT".length(), "TTTT", 'N');
    final DescriptionCommon description = TrimmingTest.getDescription(locus.getRefNts(), "TTTT:TTTC", "TTTT:TCTC");
    final VariantSample vsDad = getVariantSample(Ploidy.DIPLOID, "TTTT:TTTC", false, 5.0, VariantSample.DeNovoStatus.NOT_DE_NOVO, 0.0, description);
    final VariantSample vsMum = getVariantSample(Ploidy.DIPLOID, "TTTT:TTTC", false, 5.0, VariantSample.DeNovoStatus.NOT_DE_NOVO, 0.0, description);
    final VariantSample vsJimmy = getVariantSample(Ploidy.DIPLOID, "TTTT:TCTC", false, 5.0, VariantSample.DeNovoStatus.IS_DE_NOVO, null, description);
    final Variant v = new Variant(locus, vsDad, vsMum, vsJimmy);
    final Family fam = new Family("father", "mother", "child");
    final ChildFamilyLookup familyLookup = new ChildFamilyLookup(3, fam);
    final MendelianDenovoChecker denovoCorrector = new MendelianDenovoChecker(familyLookup);
    List<Variant> variants = Trimming.split(v, denovoCorrector);
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
    final Family fam = new Family("father", "mother", "child");
    final ChildFamilyLookup familyLookup = new ChildFamilyLookup(3, fam);
    final MendelianDenovoChecker denovoCorrector = new MendelianDenovoChecker(familyLookup);
    List<Variant> variants = Trimming.split(v, denovoCorrector);
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
    //final VariantSample vsMum = null; // new VariantSample(Ploidy.NONE, VariantSample.DeNovoStatus.UNSPECIFIED, null);
    final VariantSample vsJimmy = getVariantSample(Ploidy.HAPLOID, "TCTC", false, 5.0, VariantSample.DeNovoStatus.IS_DE_NOVO, null, description);
    final Variant v = new Variant(locus, vsDad, null, vsJimmy);
    final Family fam = new Family("father", "mother", "child");
    final ChildFamilyLookup familyLookup = new ChildFamilyLookup(3, fam);
    final MendelianDenovoChecker denovoCorrector = new MendelianDenovoChecker(familyLookup);
    List<Variant> variants = Trimming.split(v, denovoCorrector);
    final int[] expectedPositions = {1, 3};
    final String[] expectedParentA = {"T", "C"};
    final String[] expectedParentB = {null, null};
    final String[] expectedChild = {"C", "C"};
    final VariantSample.DeNovoStatus[] expectedDenovo = new VariantSample.DeNovoStatus[] {VariantSample.DeNovoStatus.IS_DE_NOVO, VariantSample.DeNovoStatus.NOT_DE_NOVO};
    checkDenovoCorrect(variants, expectedPositions, expectedParentA, expectedParentB, expectedChild, expectedDenovo);
  }

  private void checkDenovoCorrect(List<Variant> variants, int[] expectedPositions, String[] expectedParentA, String[] expectedParentB, String[] expectedChild, VariantSample.DeNovoStatus[] expectedDenovo) {
    for (int i = 0; i < variants.size(); i++) {
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
}
