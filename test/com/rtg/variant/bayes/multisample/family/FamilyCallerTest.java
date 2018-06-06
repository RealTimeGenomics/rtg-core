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

package com.rtg.variant.bayes.multisample.family;


import java.util.ArrayList;
import java.util.EnumSet;
import java.util.List;

import com.rtg.launcher.AbstractNanoTest;
import com.rtg.mode.DNA;
import com.rtg.relation.Family;
import com.rtg.relation.GenomeRelationships;
import com.rtg.relation.PedigreeException;
import com.rtg.variant.GenomePriorParams;
import com.rtg.variant.GenomePriorParamsBuilder;
import com.rtg.variant.StaticThreshold;
import com.rtg.variant.Variant;
import com.rtg.variant.Variant.VariantFilter;
import com.rtg.variant.VariantOutputLevel;
import com.rtg.variant.VariantParams;
import com.rtg.variant.VariantParamsBuilder;
import com.rtg.variant.VariantSample;
import com.rtg.variant.bayes.Description;
import com.rtg.variant.bayes.Hypotheses;
import com.rtg.variant.bayes.MockGenotypeMeasure;
import com.rtg.variant.bayes.Model;
import com.rtg.variant.bayes.ModelInterface;
import com.rtg.variant.bayes.NoAlleleBalance;
import com.rtg.variant.bayes.multisample.HaploidDiploidHypotheses;
import com.rtg.variant.bayes.multisample.HypothesisScore;
import com.rtg.variant.bayes.multisample.MultisampleJointCaller;
import com.rtg.variant.bayes.multisample.Utils;
import com.rtg.variant.bayes.snp.DescriptionSnp;
import com.rtg.variant.bayes.snp.EvidenceQ;
import com.rtg.variant.bayes.snp.HypothesesNone;
import com.rtg.variant.bayes.snp.HypothesesPrior;
import com.rtg.variant.bayes.snp.HypothesesSnp;
import com.rtg.variant.bayes.snp.StatisticsSnp;
import com.rtg.variant.format.VariantOutputVcfFormatter;
import com.rtg.variant.format.VcfFormatField;
import com.rtg.variant.util.arithmetic.SimplePossibility;

/**
 */
public class FamilyCallerTest extends AbstractNanoTest {

  private VariantOutputVcfFormatter makeFormatter(int numSamples) {
    final List<String> names = new ArrayList<>();
    for (int i = 0; i < numSamples; ++i) {
      names.add("g" + i);
    }
    return new VariantOutputVcfFormatter(names.toArray(new String[names.size()]));
  }

  protected String nanoPrefix() {
    return "familycaller-";
  }

  protected MultisampleJointCaller getFamilyCaller(Family family, VariantParams vParams) {
    return new FamilyCaller(vParams, family);
  }

  static HypothesesPrior<Description> haploidHypotheses(GenomePriorParams params, int ref) {
    return new HypothesesSnp(SimplePossibility.SINGLETON, params, true, ref - 1);
  }
  static HypothesesPrior<Description> diploidHypotheses(GenomePriorParams params, int ref) {
    return new HypothesesSnp(SimplePossibility.SINGLETON, params, false, ref - 1);
  }

  public static Family makeFamily(String father, String mother, String...children) {
    final GenomeRelationships pedigree = new GenomeRelationships();
    pedigree.addGenome(father, GenomeRelationships.SEX_MALE);
    pedigree.addGenome(mother, GenomeRelationships.SEX_FEMALE);
    for (final String child : children) {
      pedigree.addGenome(child, GenomeRelationships.SEX_MALE);
      pedigree.addParentChild(father, child);
      pedigree.addParentChild(mother, child);
    }
    //System.err.println(pedigree);
    try {
      return new Family(pedigree, father, mother, children);
    } catch (PedigreeException e) {
      throw new RuntimeException(e); // Unpossible
    }
  }

  public void testComparison() throws Exception {
    final GenomePriorParams params = getGenomePriorParams();
    final Family family = makeFamily("f", "m", "c1", "c2");
    final VariantParams vParams = new VariantParamsBuilder().maxCoverageFilter(new StaticThreshold(15)).genomePriors(params).genomeRelationships(family.pedigree()).create();
    final MultisampleJointCaller fc = getFamilyCaller(family, vParams);
    final int refNt = 1;
    final List<ModelInterface<?>> b = buildFamily(params, refNt
        , "AAAAACCCCC"
        , "AAAAA"
        , "AAGAAA"
        , "AAAAACCGCCC");
    final byte[] ref = new byte[21];
    ref[19] = 3;
    ref[20] = 1;
    final Variant v = fc.makeCall("foo", 20, 21, ref, b, new HaploidDiploidHypotheses<>(HypothesesNone.SINGLETON, haploidHypotheses(params, refNt), diploidHypotheses(params, refNt), false, null));
    final VariantOutputVcfFormatter f = new VariantOutputVcfFormatter(vParams, "f", "m", "c1", "c2");
    f.addExtraFormatFields(EnumSet.of(VcfFormatField.RQ, VcfFormatField.DN, VcfFormatField.DNP));
    assertEquals(f.formatCall(v), 13, f.formatCall(v).split("\t").length);
    mNano.check(nanoPrefix() + "comparison.vcf", f.formatCall(v), false);
  }

  private GenomePriorParams getGenomePriorParams() {
    return new GenomePriorParamsBuilder().denovoRef(0.0).denovoNonRef(0.0).contraryProbability(1).create();
  }


  static List<ModelInterface<?>> buildFamily(GenomePriorParams params, final int refNt, String... members) {
    final Hypotheses<Description> hypotheses = new HypothesesSnp(SimplePossibility.SINGLETON, params, false, refNt - 1);
    final List<ModelInterface<?>> b = new ArrayList<>();
    for (int i = 0; i < members.length; ++i) {
      b.add(new Model<>(hypotheses, new StatisticsSnp(hypotheses.description()), new NoAlleleBalance()));
      increment(b.get(i), members[i]);
    }
    for (ModelInterface<?> modelInterface : b) {
      modelInterface.freeze();
    }
    return b;
  }

  static void increment(ModelInterface<?> b, String s) {
    for (int i = 0; i < s.length(); ++i) {
      final char val = s.charAt(i);
      b.increment(new EvidenceQ(DescriptionSnp.SINGLETON, DNA.valueOf(val).ordinal() - 1, 0, 0, 0.1, 0.1, true, false, false, false));
    }
  }

  public void testComparison2() throws Exception {
    final GenomePriorParams params = getGenomePriorParams();
    final VariantParams vParams = new VariantParamsBuilder().maxCoverageFilter(new StaticThreshold(15)).genomePriors(params).create();
    final MultisampleJointCaller fc = getFamilyCaller(makeFamily("f", "m", "c1", "c2"), vParams);
    final int refNt = 1;
    final List<ModelInterface<?>> b = buildFamily(params, refNt
        , "AAAAACCCCC"
        , "CCCCCCCC"
        , "CCCCCCCCCCCCCC"
        , "AAAAACCGCCCCCCCCCCCC");
    final byte[] ref = new byte[21];
    ref[19] = 3;
    ref[20] = 1;
    final Variant v = fc.makeCall("foo", 20, 21, ref, b, new HaploidDiploidHypotheses<>(HypothesesNone.SINGLETON, haploidHypotheses(params, refNt), diploidHypotheses(params, refNt), false, null));
    final String str = makeFormatter(4).formatCall(v);
    //System.err.println(str);
    assertEquals(13, str.split("\t").length);
    mNano.check(nanoPrefix() + "comparison2", str, false);
    //assertEquals(expected, str);
  }

  public void testComparison3() throws Exception {
    final GenomePriorParams params = getGenomePriorParams();
    final VariantParams vParams = new VariantParamsBuilder().maxCoverageFilter(new StaticThreshold(15)).genomePriors(params).create();
    final MultisampleJointCaller fc = getFamilyCaller(makeFamily("f", "m", "c1", "c2"), vParams);
    final int refNt = 1;
    final List<ModelInterface<?>> b = buildFamily(params, refNt
        , "AAAAAAA"
        , "CCCCCCAAA"
        , "AAAAAAA"
        , "AAAAAAA");
    final byte[] ref = new byte[21];
    ref[19] = 3;
    ref[20] = 1;
    final Variant v = fc.makeCall("foo", 20, 21, ref, b, new HaploidDiploidHypotheses<>(HypothesesNone.SINGLETON, haploidHypotheses(params, refNt), diploidHypotheses(params, refNt), false, null));
    assertEquals(13, makeFormatter(4).formatCall(v).split("\t").length);
    mNano.check(nanoPrefix() + "comparison3", makeFormatter(4).formatCall(v), false);
  }

  public void testComparison4() throws Exception {
    final GenomePriorParams params = getGenomePriorParams();
    final VariantParams vParams = new VariantParamsBuilder().maxCoverageFilter(new StaticThreshold(15)).genomePriors(params).create();
    final MultisampleJointCaller fc = getFamilyCaller(makeFamily("f", "m", "c1", "c2"), vParams);
    final int refNt = 1;
    final List<ModelInterface<?>> b = buildFamily(params, refNt
        , "CCCCCCCAAAA"
        , "AAAAAAA"
        , "AAAAAAA"
        , "AAAAAAAAA");
    final byte[] ref = new byte[21];
    ref[19] = 3;
    ref[20] = 1;
    final Variant v = fc.makeCall("foo", 20, 21, ref, b, new HaploidDiploidHypotheses<>(HypothesesNone.SINGLETON, haploidHypotheses(params, refNt), diploidHypotheses(params, refNt), false, null));
    assertEquals(13, makeFormatter(4).formatCall(v).split("\t").length);
    mNano.check(nanoPrefix() + "comparison4", makeFormatter(4).formatCall(v), false);
  }

  public void testAllOutput() throws Exception {
    final GenomePriorParams params = getGenomePriorParams();
    final VariantParams vParams = new VariantParamsBuilder()
    .callLevel(VariantOutputLevel.ALL)
    .maxCoverageFilter(new StaticThreshold(15))
    .genomePriors(params)
    .create();
    final MultisampleJointCaller fc = getFamilyCaller(makeFamily("f", "m", "c1", "c2"), vParams);
    final int refNt = 1;
    final List<ModelInterface<?>> b = buildFamily(params, refNt
        , "AAAAAAA"
        , "AAAAAAA"
        , "AAAAAAA"
        , "AAAAAAA");
    final byte[] ref = new byte[21];
    ref[19] = 3;
    ref[20] = 1;
    final Variant v = fc.makeCall("foo", 20, 21, ref, b, new HaploidDiploidHypotheses<>(HypothesesNone.SINGLETON, haploidHypotheses(params, refNt), diploidHypotheses(params, refNt), false, null));
    mNano.check(nanoPrefix() + "allOutput", makeFormatter(4).formatCall(v), false);
  }

  // Test where all family members agree, but differ from reference -- was a bug in old short-circuit code
  public void testComparison5() throws Exception {
    final GenomePriorParams params = getGenomePriorParams();
    final VariantParams vParams = new VariantParamsBuilder().maxCoverageFilter(new StaticThreshold(15)).genomePriors(params).create();
    final MultisampleJointCaller fc = getFamilyCaller(makeFamily("f", "m", "c1", "c2"), vParams);
    final int refNt = 1;
    final List<ModelInterface<?>> b = buildFamily(params, refNt
        , "CCCCCCCCC"
        , "CCCCCCCCC"
        , "CCCCCCC"
        , "CCCCCCC");
    final byte[] ref = new byte[21];
    ref[19] = 3; //G previous nt
    ref[20] = 1; //A
    final Variant v = fc.makeCall("foo", 20, 21, ref, b, new HaploidDiploidHypotheses<>(HypothesesNone.SINGLETON, haploidHypotheses(params, refNt), diploidHypotheses(params, refNt), false, null));
    assertEquals(13, makeFormatter(4).formatCall(v).split("\t").length);
    mNano.check(nanoPrefix() + "comparison5", makeFormatter(4).formatCall(v), false);
  }

  // Test where all family members agree, but differ from reference -- was a bug in old short-circuit code
  public void testComparison6() throws Exception {
    final GenomePriorParams params = getGenomePriorParams();
    final VariantParams vParams = new VariantParamsBuilder().maxCoverageFilter(new StaticThreshold(100)).genomePriors(params).create();
    final MultisampleJointCaller fc = getFamilyCaller(makeFamily("f", "m", "c1"), vParams);
    final int refNt = 4;
    final List<ModelInterface<?>> b = buildFamily(params, refNt
        , "GGGGGGGGGG" + "GGGGGGGGGG" + "GGGGGGGGGG" + "A"
            , "GGGGGGGGGG" + "GGGGGGGGGG" + "G"
                , "GGGGGGGGGG" + "GGGGGGGGGG" + "GGGGGGGGGG" + "G");
    final byte[] ref = new byte[21];
    ref[19] = 1; //A previous nt
    ref[20] = 4; //T
    final Variant v = fc.makeCall("foo", 20, 21, ref, b, new HaploidDiploidHypotheses<>(HypothesesNone.SINGLETON, haploidHypotheses(params, refNt), diploidHypotheses(params, refNt), false, null));
    final String str = makeFormatter(3).formatCall(v);
    //System.err.println(str);
    assertEquals(12, str.split("\t").length);
    assertFalse(v.isFiltered(VariantFilter.COVERAGE));
    assertEquals('A', v.getLocus().getPreviousRefNt());
    mNano.check(nanoPrefix() + "comparison6", makeFormatter(3).formatCall(v), false);
  }

  public void testOverCoverage() {
    final GenomePriorParams params = getGenomePriorParams();
    final VariantParams vParams = new VariantParamsBuilder().maxCoverageFilter(new StaticThreshold(5)).genomePriors(params).create();
    final MultisampleJointCaller fc = getFamilyCaller(makeFamily("f", "m", "c1"), vParams);
    final int refNt = 4;
    final List<ModelInterface<?>> b = buildFamily(params, refNt
        , "GGGGGGGGGG" + "GGGGGGGGGG" + "GGGGGGGGGG" + "A"
            , "GGGGGGGGGG" + "GGGGGGGGGG" + "G"
                , "GGGGGGGGGG" + "GGGGGGGGGG" + "GGGGGGGGGG" + "G");
    final byte[] ref = new byte[21];
    ref[19] = 1; //A previous nt
    ref[20] = 4; //T
    final Variant v = fc.makeCall("foo", 20, 21, ref, b, new HaploidDiploidHypotheses<>(HypothesesNone.SINGLETON, haploidHypotheses(params, refNt), diploidHypotheses(params, refNt), false, null));
    assertTrue(v.isFiltered(VariantFilter.COVERAGE));
  }

  public void testHardLeft() {
    final GenomePriorParams params = getGenomePriorParams();
    final VariantParams vParams = new VariantParamsBuilder().maxCoverageFilter(new StaticThreshold(5)).genomePriors(params).create();
    final MultisampleJointCaller fc = getFamilyCaller(makeFamily("f", "m", "c1"), vParams);
    final int refNt = 4;
    final List<ModelInterface<?>> b = buildFamily(params, refNt
        , "GGGGGGGGGG" + "GGGGGGGGGG" + "GGGGGGGGGG" + "A"
            , "GGGGGGGGGG" + "GGGGGGGGGG" + "G"
                , "GGGGGGGGGG" + "GGGGGGGGGG" + "GGGGGGGGGG" + "G");
    final byte[] ref = new byte[21];
    ref[0] = 4; //T previous nt
    final Variant v = fc.makeCall("foo", 0, 1, ref, b, new HaploidDiploidHypotheses<>(HypothesesNone.SINGLETON, haploidHypotheses(params, refNt), diploidHypotheses(params, refNt), false, null));
    assertEquals('N', v.getLocus().getPreviousRefNt());
  }

  public void testMeanCoverage() {
    final GenomePriorParams params = getGenomePriorParams();
    final int refNt = 4;
    final List<ModelInterface<?>> models = buildFamily(params, refNt , "GGGGGGGGGG" + "GGGGGTTTTT" + "GGGGGGGGGG" + "A", "GGGGGGGGGG" + "GGGGGGGGGG" + "G", "GGGCCCCGGG" + "GGGGGGGGGG" + "GGGGGGGGGG" + "G");
    assertEquals(31, Utils.maxCoverage(models));
    assertEquals(0, Utils.maxCoverage(new ArrayList<>()));
    assertEquals(0.0, Utils.meanAmbiguityRatio(models));
  }

  public void testComparisonPloidy() throws Exception {
    final GenomePriorParams params = getGenomePriorParams();
    final VariantParams vParams = new VariantParamsBuilder().maxCoverageFilter(new StaticThreshold(15)).genomePriors(params).create();
    final int refNt = 1;
    final List<ModelInterface<?>> b = buildFamily(params, refNt
        , "CCCCCCCAAAA"
        , "AAAAAAA"
        , "AAAAAAA"
        , "AAAAAAAAA");
    final Hypotheses<Description> hypotheses = haploidHypotheses(params, 1);
    Model<Description> model = new Model<>(hypotheses, new StatisticsSnp(hypotheses.description()), new NoAlleleBalance());
    increment(model, "CCCCC");
    model.freeze();
    b.add(Family.FATHER_INDEX, model);
    model = new Model<>(hypotheses, new StatisticsSnp(hypotheses.description()), new NoAlleleBalance());
    increment(model, "AACCCCC");
    model.freeze();
    b.add(model);
    final MultisampleJointCaller fc = getFamilyCaller(makeFamily("f", "m", "c1", "c2", "c3", "c4"), vParams);
    final byte[] ref = new byte[21];
    ref[19] = 3;
    ref[20] = 1;
    final Variant v = fc.makeCall("foo", 20, 21, ref, b, new HaploidDiploidHypotheses<>(HypothesesNone.SINGLETON, haploidHypotheses(params, refNt), diploidHypotheses(params, refNt), false, null));
    assertEquals(15, makeFormatter(6).formatCall(v).split("\t").length);
    mNano.check(nanoPrefix() + "comparisonPloidy", makeFormatter(6).formatCall(v), false);
  }
  public void testBoring() {
    final GenomePriorParams params = getGenomePriorParams();
    final VariantParams vParams = new VariantParamsBuilder().maxCoverageFilter(new StaticThreshold(15)).genomePriors(params).create();
    final MultisampleJointCaller fc = getFamilyCaller(makeFamily("f", "m", "c1", "c2"), vParams);
    final int refNt = 1;
    final List<ModelInterface<?>> b = buildFamily(params, refNt
        , "CCAAAAAAAAA"
        , "AAAAAAA"
        , "AAAAAAA"
        , "AAAAAAAAA");
    final byte[] ref = new byte[21];
    ref[19] = 3;
    ref[20] = 1;
    final Variant v = fc.makeCall("foo", 20, 21, ref, b, new HaploidDiploidHypotheses<>(HypothesesNone.SINGLETON, haploidHypotheses(params, refNt), diploidHypotheses(params, refNt), false, null));
    assertNull(v);
  }

  public void testDisagreement() {
    final HypothesisScore newScore = new HypothesisScore(new MockGenotypeMeasure(0, 0, 5, 0));

    final HypothesisScore[] scores = new HypothesisScore[3];
    scores[0] = new HypothesisScore(new MockGenotypeMeasure(0, 1, 10, 0));
    scores[1] = new HypothesisScore(new MockGenotypeMeasure(0, 1, 2, 0));
    scores[2] = new HypothesisScore(new MockGenotypeMeasure(0, 0, 7, 0));
    for (HypothesisScore score : scores) {
      score.setDenovo(VariantSample.DeNovoStatus.IS_DE_NOVO);
      score.setDeNovoPosterior(0.3);
    }
    assertFalse(FamilyCaller.setScore(scores, 0, newScore));
    assertFalse(FamilyCaller.setScore(scores, 1, newScore));
    assertTrue(FamilyCaller.setScore(scores, 2, newScore));

  }

}

