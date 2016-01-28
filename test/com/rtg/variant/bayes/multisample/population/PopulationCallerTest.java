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

package com.rtg.variant.bayes.multisample.population;


import static com.rtg.util.StringUtils.LS;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.StringReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.rtg.mode.DNA;
import com.rtg.relation.Family;
import com.rtg.relation.GenomeRelationships;
import com.rtg.relation.MultiFamilyOrdering;
import com.rtg.relation.PedigreeException;
import com.rtg.util.InvalidParamsException;
import com.rtg.variant.GenomePriorParams;
import com.rtg.variant.GenomePriorParamsBuilder;
import com.rtg.variant.StaticThreshold;
import com.rtg.variant.Variant;
import com.rtg.variant.VariantOutputLevel;
import com.rtg.variant.VariantParams;
import com.rtg.variant.VariantParamsBuilder;
import com.rtg.variant.VariantSample;
import com.rtg.variant.bayes.Description;
import com.rtg.variant.bayes.Model;
import com.rtg.variant.bayes.ModelInterface;
import com.rtg.variant.bayes.ModelNone;
import com.rtg.variant.bayes.NoAlleleBalance;
import com.rtg.variant.bayes.multisample.HaploidDiploidHypotheses;
import com.rtg.variant.bayes.multisample.HypothesisScore;
import com.rtg.variant.bayes.multisample.HypothesisScores;
import com.rtg.variant.bayes.multisample.PriorContainer;
import com.rtg.variant.bayes.multisample.Utils;
import com.rtg.variant.bayes.multisample.family.AbstractFamilyPosteriorTest;
import com.rtg.variant.bayes.multisample.family.FamilyCaller;
import com.rtg.variant.bayes.multisample.forwardbackward.FamilyCallerFB;
import com.rtg.variant.bayes.multisample.forwardbackward.FamilyCallerFBTest;
import com.rtg.variant.bayes.snp.DescriptionCommon;
import com.rtg.variant.bayes.snp.DescriptionSnp;
import com.rtg.variant.bayes.snp.EvidenceQ;
import com.rtg.variant.bayes.snp.HypothesesNone;
import com.rtg.variant.bayes.snp.HypothesesPrior;
import com.rtg.variant.bayes.snp.HypothesesSnp;
import com.rtg.variant.bayes.snp.StatisticsSnp;
import com.rtg.variant.format.VariantOutputVcfFormatter;
import com.rtg.variant.util.arithmetic.LogApproximatePossibility;
import com.rtg.variant.util.arithmetic.LogPossibility;
import com.rtg.variant.util.arithmetic.PossibilityArithmetic;

import junit.framework.TestCase;

/**
 *         Date: 8/03/12
 *         Time: 11:36 AM
 */
public class PopulationCallerTest extends TestCase {

  static final int PRIOR_FIELDS = 9;
  static final int FORMAT = 8;

  public static HypothesesSnp haploidHypotheses(GenomePriorParams params, int refNt) {
    return new HypothesesSnp(LogPossibility.SINGLETON, params, true, refNt - 1);
  }

  public static HypothesesSnp diploidHypotheses(GenomePriorParams params, int refNt) {
    return new HypothesesSnp(LogPossibility.SINGLETON, params, false, refNt - 1);
  }

  public static List<ModelInterface<?>> buildModels(GenomePriorParams params, final int refNt, String... members) {
    final HypothesesSnp hypotheses = new HypothesesSnp(LogPossibility.SINGLETON, params, false, refNt - 1);
    return buildModels(hypotheses, members);
  }

  private static List<ModelInterface<?>> buildModels(HypothesesSnp hypotheses, String... members) {
    final List<ModelInterface<?>> b = new ArrayList<>();
    for (final String member : members) {
      final Model<Description> model = new Model<>(hypotheses, new StatisticsSnp(hypotheses.description()), new NoAlleleBalance());
      b.add(model);
      increment(model, member, 0.1);
    }
    return b;
  }

  public static <T extends Description> void increment(ModelInterface<T> b, String s, double r) {
    increment(b, s, r, 0.1);
  }
  public static <T extends Description> void increment(ModelInterface<T> b, String s, double r, double q) {
    for (int i = 0; i < s.length(); i++) {
      final char val = s.charAt(i);
      b.increment(new EvidenceQ(DescriptionSnp.SINGLETON, DNA.valueOf(val).ordinal() - 1, 0, 0, r, q, true, false, false, false));
    }
  }

  public static VariantOutputVcfFormatter makeFormatter(VariantParams params, int numSamples) {
    final List<String> names = new ArrayList<>();
    for (int i = 0; i < numSamples; i++) {
      names.add("g" + i);
    }
    final String[] namesarr = names.toArray(new String[names.size()]);
    if (params != null) {
      return new VariantOutputVcfFormatter(params, namesarr);
    }
    return new VariantOutputVcfFormatter(namesarr);
  }

  public void testComparison() throws InvalidParamsException, IOException {
    final GenomePriorParams params = new GenomePriorParamsBuilder().create();
    final VariantParams vParams = new VariantParamsBuilder().maxCoverageFilter(new StaticThreshold(15)).genomeRelationships(new GenomeRelationships()).create();
    final PopulationCaller pc = new PopulationCaller(vParams);
    final int refNt = 1;
    final List<ModelInterface<?>> b = buildModels(params, refNt
        , "AAAAACCCCC"
        , "CCCCCCCC"
        , "CCCCCCCCCCCCCC"
        , "AAAAACCGCCCCCCCCCCCC");
    final byte[] ref = new byte[21];
    ref[19] = 3;
    ref[20] = 1;
    final Variant v = pc.makeCall("foo", 20, 21, ref, b, new HaploidDiploidHypotheses<>(HypothesesNone.SINGLETON, haploidHypotheses(params, refNt), diploidHypotheses(params, refNt)));
    final String str = makeFormatter(vParams, 4).formatCall(v);
    //System.err.println(str);
    assertEquals(13, str.split("\t").length);
    final String[] parts = str.split("\t");
    assertEquals("21", parts[1]);
    assertEquals("OC", parts[6]);
    // Following numbers not actually calculated externally
    testCall(str, "GQ", "30", "23", "40", "7");
  }

  public void testComparisonAmbiguousRatio() throws InvalidParamsException, IOException {
    final GenomePriorParams params = new GenomePriorParamsBuilder().create();
    final VariantParams vParams = new VariantParamsBuilder().maxAmbiguity(0.03).maxCoverageFilter(new StaticThreshold(50)).genomePriors(params).genomeRelationships(new GenomeRelationships()).create();
    final PopulationCaller pc = new PopulationCaller(vParams);
    final int refNt = 1;
    final List<ModelInterface<?>> b = buildModels(params, refNt
        , "AAAAACCCCC"
        , "CCCCCCCC"
        , "CCCCCCCCCCCCCC"
        , "AAAAACCGCCCCCCCCCCCC");
    increment(b.get(0), "A", 0.6);
    final byte[] ref = new byte[21];
    ref[19] = 3;
    ref[20] = 1;
    final Variant v = pc.makeCall("foo", 20, 21, ref, b, new HaploidDiploidHypotheses<>(HypothesesNone.SINGLETON, haploidHypotheses(params, refNt), diploidHypotheses(params, refNt)));
    final String str = makeFormatter(vParams, 4).formatCall(v);
    //System.err.println(str);
    assertEquals(13, str.split("\t").length);
    final String[] parts = str.split("\t");
    assertEquals("PASS", parts[6]);
    testCall(str, "AR", "0.091", "0.000", "0.000", "0.000");
  }

  public void testN() throws InvalidParamsException, IOException {
    final GenomePriorParams params = new GenomePriorParamsBuilder().create();
    final VariantParams vParams = new VariantParamsBuilder().maxCoverageFilter(new StaticThreshold(15)).genomePriors(params).create();
    final PopulationCaller pc = new PopulationCaller(vParams);
    final int refNt = 1;
    final List<ModelInterface<?>> b = buildModels(params, refNt
        , "GGGGG"
        , "CCCCCCCCC"
        , "GGGGGAAAAAAAAA"
        , "CCCCCAAAAAAAAAAAAAAA");
    final byte[] ref = new byte[21];
    ref[19] = 3;
    ref[20] = 0;
    final Variant v = pc.makeCall("foo", 20, 21, ref, b, new HaploidDiploidHypotheses<>(HypothesesNone.SINGLETON, haploidHypotheses(params, refNt), diploidHypotheses(params, refNt)));
    assertNotNull(v);
  }

  public void testShortCircuit() throws InvalidParamsException, IOException {
    final GenomePriorParams params = new GenomePriorParamsBuilder().create();
    final VariantParams vParams = new VariantParamsBuilder().maxCoverageFilter(new StaticThreshold(15)).genomePriors(params).create();
    final PopulationCaller pc = new PopulationCaller(vParams);
    final int refNt = 1;
    final List<ModelInterface<?>> b = buildModels(params, refNt
        , "AAAAAAAAAA"
        , "AAAAAAAA"
        , "AAAAAAAAAAAAAA"
        , "AAAAAAAAAAAAAAAAAAAA");
    final byte[] ref = new byte[21];
    ref[19] = 3;
    ref[20] = 1;
    final Variant v = pc.makeCall("foo", 20, 21, ref, b, new HaploidDiploidHypotheses<>(HypothesesNone.SINGLETON, haploidHypotheses(params, refNt), diploidHypotheses(params, refNt)));
    assertNull(v);
  }

  public void testUninteresting() throws InvalidParamsException, IOException {
    final GenomePriorParams params = new GenomePriorParamsBuilder().create();
    final VariantParams vParams = new VariantParamsBuilder().maxCoverageFilter(new StaticThreshold(15)).genomePriors(params).create();
    final PopulationCaller pc = new PopulationCaller(vParams);
    final int refNt = 1;
    final List<ModelInterface<?>> b = buildModels(params, refNt
        , "AAAAAAAAAA"
        , "AAAAAAAC"
        , "AAAAAAACAAAAAA"
        , "AAAAAAAAACAAAAACAAAA");
    final byte[] ref = new byte[21];
    ref[19] = 3;
    ref[20] = 1;
    final Variant v = pc.makeCall("foo", 20, 21, ref, b, new HaploidDiploidHypotheses<>(HypothesesNone.SINGLETON, haploidHypotheses(params, refNt), diploidHypotheses(params, refNt)));
    assertNull(v);
  }
  public void testAllShortCurcuit() throws InvalidParamsException, IOException {
    final GenomePriorParams params = new GenomePriorParamsBuilder().create();
    final VariantParams vParams = new VariantParamsBuilder().callLevel(VariantOutputLevel.ALL).maxCoverageFilter(new StaticThreshold(15)).genomePriors(params).create();
    final PopulationCaller pc = new PopulationCaller(vParams);
    final int refNt = 1;
    final List<ModelInterface<?>> b = buildModels(params, refNt
        , "AAAAAAAAAA"
        , "AAAAAAAA"
        , "AAAAAAAAAAAAAA"
        , "AAAAAAAAAAAAAAAAAAAA");
    final byte[] ref = new byte[21];
    ref[19] = 3;
    ref[20] = 1;
    final Variant v = pc.makeCall("foo", 20, 21, ref, b, new HaploidDiploidHypotheses<>(HypothesesNone.SINGLETON, haploidHypotheses(params, refNt), diploidHypotheses(params, refNt)));
    assertNotNull(v);
    assertFalse(v.isInteresting());
  }

  public void testAll() throws InvalidParamsException, IOException {
    final GenomePriorParams params = new GenomePriorParamsBuilder().create();
    final VariantParams vParams = new VariantParamsBuilder().callLevel(VariantOutputLevel.ALL).maxCoverageFilter(new StaticThreshold(15)).genomePriors(params).create();
    final PopulationCaller pc = new PopulationCaller(vParams);
    final int refNt = 1;
    final List<ModelInterface<?>> b = buildModels(params, refNt
        , "AAAAAAAAAA"
        , "AAAAAAAA"
        , "AAAAAAAAAAAAAA"
        , "AAAAAAAAAAAAAAGGGGGG");
    final byte[] ref = new byte[21];
    ref[19] = 3;
    ref[20] = 1;
    final Variant v = pc.makeCall("foo", 20, 21, ref, b, new HaploidDiploidHypotheses<>(HypothesesNone.SINGLETON, haploidHypotheses(params, refNt), diploidHypotheses(params, refNt)));
    assertNotNull(v);
    assertFalse(v.isInteresting());
    assertNotNull(v.getNonIdentityPosterior());
    assertTrue(v.getNonIdentityPosterior() < 0.0);
  }

  public void testAllInteresting() throws InvalidParamsException, IOException {
    final GenomePriorParams params = new GenomePriorParamsBuilder().create();
    final VariantParams vParams = new VariantParamsBuilder().callLevel(VariantOutputLevel.ALL).maxCoverageFilter(new StaticThreshold(15)).genomePriors(params).create();
    final PopulationCaller pc = new PopulationCaller(vParams);
    final int refNt = 1;
    final List<ModelInterface<?>> b = buildModels(params, refNt
        , "AAAAAAAAAA"
        , "AAAAAAAA"
        , "AAAAAAAAAAAAAA"
        , "GGGGGGGGGGGGAAGGGGGG");
    final byte[] ref = new byte[21];
    ref[19] = 3;
    ref[20] = 1;
    final Variant v = pc.makeCall("foo", 20, 21, ref, b, new HaploidDiploidHypotheses<>(HypothesesNone.SINGLETON, haploidHypotheses(params, refNt), diploidHypotheses(params, refNt)));
    assertNotNull(null, v);
    assertTrue(v.isInteresting());
    assertNotNull(v.getNonIdentityPosterior());
    assertTrue(v.getNonIdentityPosterior() > 0.0);
  }

  public void testAmbiguity() throws InvalidParamsException, IOException {
    final GenomePriorParams params = new GenomePriorParamsBuilder().create();
    final VariantParams vParams = new VariantParamsBuilder().callLevel(VariantOutputLevel.ALL).maxAmbiguity(0.3).maxCoverageFilter(new StaticThreshold(50)).genomePriors(params).genomeRelationships(new GenomeRelationships()).create();
    final PopulationCaller pc = new PopulationCaller(vParams);
    final int refNt = 1;
    final List<ModelInterface<?>> b = buildModels(params, refNt
        , "AAAAGGGGGGG"
        , "AAAAAAAA");
    final HypothesesSnp hypotheses = new HypothesesSnp(LogPossibility.SINGLETON, params, false, refNt - 1);
    Model<Description> model = new Model<>(hypotheses, new StatisticsSnp(hypotheses.description()), new NoAlleleBalance());
    increment(model, "GGGCGG", 0.6);
    increment(model, "GGGCGG", 0.1);
    b.add(model);
    model = new Model<>(hypotheses, new StatisticsSnp(hypotheses.description()), new NoAlleleBalance());
    increment(model, "GGCGGCCCGGGG", 0.6);
    increment(model, "GGGC", 0.1);
    b.add(model);
    final byte[] ref = new byte[21];
    ref[19] = 3;
    ref[20] = 1;
    final Variant v = pc.makeCall("foo", 20, 21, ref, b, new HaploidDiploidHypotheses<>(HypothesesNone.SINGLETON, haploidHypotheses(params, refNt), diploidHypotheses(params, refNt)));
    assertNotNull(null, v);
    assertTrue(v.isInteresting());
    final String str = makeFormatter(vParams, 4).formatCall(v);
    assertEquals(13, str.split("\t").length);
    assertTrue(str.contains("a30.0"));

    assertNotNull(v.getNonIdentityPosterior());
  }

  public void testBaseCoverage() {
    assertEquals(0, Utils.maxCoverage(new ArrayList<ModelInterface<?>>()));
  }

  public void testPloidy() throws InvalidParamsException, IOException {
    final GenomePriorParams params = new GenomePriorParamsBuilder().create();
    final VariantParams vParams = new VariantParamsBuilder().callLevel(VariantOutputLevel.ALL).maxAmbiguity(0.3).maxCoverageFilter(new StaticThreshold(50)).genomePriors(params).genomeRelationships(new GenomeRelationships()).create();
    final PopulationCaller pc = new PopulationCaller(vParams);
    final int refNt = 1;
    final List<ModelInterface<?>> b = buildModels(params, refNt
        , "AAAAGGGGGGG"
        , "AAAAAAAA");
    final HypothesesSnp hypotheses = new HypothesesSnp(LogPossibility.SINGLETON, params, true, refNt - 1);
    Model<Description> model = new Model<>(hypotheses, new StatisticsSnp(hypotheses.description()), new NoAlleleBalance());
    increment(model, "GGGCGG", 0.6);
    increment(model, "GGGCGG", 0.1);
    b.add(model);
    model = new Model<>(hypotheses, new StatisticsSnp(hypotheses.description()), new NoAlleleBalance());
    increment(model, "GGCGGCCCGGGG", 0.6);
    increment(model, "GGGC", 0.1);
    b.add(model);
    final byte[] ref = new byte[21];
    ref[19] = 3;
    ref[20] = 1;
    final Variant v = pc.makeCall("foo", 20, 21, ref, b, new HaploidDiploidHypotheses<>(HypothesesNone.SINGLETON, haploidHypotheses(params, refNt), diploidHypotheses(params, refNt)));
    assertNotNull(null, v);
    assertTrue(v.isInteresting());
    final String str = makeFormatter(vParams, 4).formatCall(v);
    assertEquals(13, str.split("\t").length);
    testCall(str,  "GT", "0/1", "0/0", "1", "1");
  }


  public static int findField(String output, String name) {
    final String[] format = output.split("\t")[FORMAT].split(":");
    return Arrays.asList(format).indexOf(name);
  }


  public static void testCall(String output, String fieldName, String... calls) {
    final int fieldIndex = findField(output, fieldName);
    final String[] parts = output.trim().split("\t");
    assertEquals(calls.length, parts.length - PRIOR_FIELDS);
    for (int i = 0; i < parts.length - PRIOR_FIELDS; i++) {
      final String[] subParts = parts[i + PRIOR_FIELDS].split(":");
      assertEquals("The field \"" + fieldName + " was mismatched in: <" + output + "> expected values <" + Arrays.toString(calls) + ">", calls[i], subParts[fieldIndex]);
    }
  }

  public void testDisagreeingRefInteresting() {
    // See bug 1534. The root cause was that a family call within a population could be different from reference but
    // be overridden by a call for another family. At the time of the bug the "interesting" property was being
    // recomputed within the population caller by checking for calls in the final result that didn't match the reference
    // but this recomputed value was not used consistently. An adjacent complex region was not extended to cover this
    // call but then had a previous base prepended resulting in an overlap.
    final GenomePriorParams params = new GenomePriorParamsBuilder().contraryProbability(1).create();
    final VariantParams vParams = new VariantParamsBuilder().callLevel(VariantOutputLevel.ALL).maxAmbiguity(0.3).maxCoverageFilter(new StaticThreshold(50)).genomePriors(params).create();
    final Family[] families = new Family[2];
    families[0] = new Family("0", "1", "2");
    families[0].setSampleId(0, 0);
    families[0].setSampleId(1, 1);
    families[0].setSampleId(2, 2);
    families[1] = new Family("2", "3", "4", "5", "6");
    families[1].setSampleId(0, 2);
    families[1].setSampleId(1, 3);
    families[1].setSampleId(2, 4);
    families[1].setSampleId(3, 5);
    families[1].setSampleId(4, 6);
    final FamilyCaller familyCaller = new FamilyCaller(vParams, families);
    final PopulationCaller pc = new PopulationCaller(vParams, familyCaller);
    final int refNt = 2;
    final List<ModelInterface<?>> b = new ArrayList<>();
    final HypothesesSnp hypothesesSnp =  new HypothesesSnp(LogApproximatePossibility.SINGLETON, params, true, refNt - 1);
    final HypothesesSnp diploid =  new HypothesesSnp(LogApproximatePossibility.SINGLETON, params, false, refNt - 1);

    final String[] values = {
        "CCCCCCCCCCCCGGGGGG"
        , ""
        , "CCCCGGGGGGGG"
        , ""
        , "CCCCCCCC"
        , "CCCCCCCC"
        , "CCCCCCCC"
    };
    for (int i = 0; i < 7; i++) {
      final Model<Description> model = new Model<>(hypothesesSnp, new StatisticsSnp(hypothesesSnp.description()), new NoAlleleBalance());
      b.add(model);
      increment(model, values[i], 0.001, 0.001);
    }
    // Make it haploid
    b.set(1, ModelNone.SINGLETON);
    b.set(3, ModelNone.SINGLETON);
    final byte[] ref = new byte[21];
    ref[19] = 3;
    ref[20] = 1;
    final Variant v = pc.makeCall("foo", 20, 21, ref, b, new HaploidDiploidHypotheses<>(HypothesesNone.SINGLETON, hypothesesSnp, diploid));
    assertTrue(v.toString(), v.isInteresting());
    /* if USE_FB_FALLBACK
    for (int i = 0; i < v.getNumberOfSamples(); i++) {
      final VariantSample sample = v.getSample(i);
      if (sample != null) {
        assertEquals("C", sample.getName());
      }
    }
    // No longer interesting now that it is deferred to the FB caller
    assertFalse(v.toString(), v.isInteresting());
   */
  }

  /*
      Fathers left, Mothers right


                gdad---gmom
                     |
        dad---------mom-----dad2
            |  |  |     |  |
           S1 D1 S2     S3 D2
      Attempt to test pedigree calling on Y chromosome
   */
  public void testPedigreeHaploidNoneHaploid2() throws InvalidParamsException {
    final GenomePriorParams priors = new GenomePriorParamsBuilder().denovoRef(0.01).denovoNonRef(0.00001).create();

    final GenomeRelationships pedigree = new GenomeRelationships();
    pedigree.addGenome("gdad", GenomeRelationships.SEX_MALE);
    pedigree.addGenome("gmom", GenomeRelationships.SEX_FEMALE);
    pedigree.addGenome("dad", GenomeRelationships.SEX_MALE);
    pedigree.addGenome("mom", GenomeRelationships.SEX_FEMALE);
    pedigree.addGenome("dad2", GenomeRelationships.SEX_MALE);
    pedigree.addGenome("S1", GenomeRelationships.SEX_MALE);
    pedigree.addGenome("D1", GenomeRelationships.SEX_FEMALE);
    pedigree.addGenome("S2", GenomeRelationships.SEX_MALE);
    pedigree.addGenome("S3", GenomeRelationships.SEX_MALE);
    pedigree.addGenome("D2", GenomeRelationships.SEX_FEMALE);

    pedigree.addParentChild("gdad", "mom");
    pedigree.addParentChild("gmom", "mom");
    pedigree.addParentChild("dad", "S1");
    pedigree.addParentChild("dad", "D1");
    pedigree.addParentChild("dad", "S2");
    pedigree.addParentChild("mom", "S1");
    pedigree.addParentChild("mom", "D1");
    pedigree.addParentChild("mom", "S2");
    pedigree.addParentChild("mom", "S3");
    pedigree.addParentChild("mom", "D2");
    pedigree.addParentChild("dad2", "S3");
    pedigree.addParentChild("dad2", "D2");

    final PossibilityArithmetic arith = LogApproximatePossibility.SINGLETON;
    final Description desc = new DescriptionCommon("A", "C");

    final HypothesesPrior<Description> none = HypothesesNone.SINGLETON;
    final HypothesesPrior<Description> haploid = new AbstractFamilyPosteriorTest.MockHyp(desc, arith, true, new double[] {0.2, 0.8});
    final HaploidDiploidHypotheses<HypothesesPrior<Description>> hdh = new HaploidDiploidHypotheses<>(none, haploid, null, false, null);

    final List<ModelInterface<?>> models = new ArrayList<>();
    models.add(new AbstractFamilyPosteriorTest.MockModel(haploid, new double[]{0.25, 0.75}));
    models.add(new AbstractFamilyPosteriorTest.MockModel(none));
    models.add(new AbstractFamilyPosteriorTest.MockModel(haploid, new double[]{0.3, 0.7}));
    models.add(new AbstractFamilyPosteriorTest.MockModel(none));
    models.add(new AbstractFamilyPosteriorTest.MockModel(haploid, new double[]{0.3, 0.7}));
    models.add(new AbstractFamilyPosteriorTest.MockModel(haploid, new double[]{0.3, 0.7}));
    models.add(new AbstractFamilyPosteriorTest.MockModel(none));
    models.add(new AbstractFamilyPosteriorTest.MockModel(haploid, new double[]{0.25, 0.75}));
    models.add(new AbstractFamilyPosteriorTest.MockModel(haploid, new double[]{0.25, 0.75}));
    models.add(new AbstractFamilyPosteriorTest.MockModel(none));

    final Set<Family> families = Family.getFamilies(pedigree, false, null);
    final List<String> calledGenomes = new ArrayList<>();
    // Ensure all family members are in called genomes
    for (Family family : families) {
      for (String member : family.getMembers()) {
        if (!calledGenomes.contains(member)) {
          calledGenomes.add(member);
        }
      }
    }
    for (Family family : families) {
      family.setSampleIds(calledGenomes);
    }
    final FamilyCallerFB ffb = new FamilyCallerFB(VariantParams.builder().maxEmIterations(0).genomePriors(priors).create()
        , families.toArray(new Family[families.size()])
    );
    final PopulationCaller popCall = new PopulationCaller(VariantParams.builder().genomePriors(priors).create(), ffb);
    final HypothesisScores bestScores = popCall.getBestScores(models, new PriorContainer<>(hdh, ffb.makeInitialBs(models)));
    assertFalse(Double.isNaN(bestScores.getNonIdentityPosterior()));
    assertFalse(Double.isInfinite(bestScores.getNonIdentityPosterior()));

    final HypothesisScore[] scores = bestScores.getScores();
    FamilyCallerFBTest.checkScores(scores);

  }

  public void testDisagreeingPPPFallbackEM() {
    // Brian saw a crash when running with EM iterations. The family caller had to fall back to PPP which caused a null pointer exception in FB when accessing the updated priors.
    // FamilyCaller was returning a HypothesisScore with null Bs
    // This happens when the initial calls agree but the second EM iteration results in a disagreement
    final GenomePriorParams params = new GenomePriorParamsBuilder().create();
    final VariantParams vParams = new VariantParamsBuilder()
        .callLevel(VariantOutputLevel.ALL)
        .maxAmbiguity(0.3)
        .maxCoverageFilter(new StaticThreshold(50))
        .genomePriors(params)
        .maxEmIterations(6)
        .create();
    final Family[] families = new Family[3];
    families[0] = new Family("0", "1", "2");
    families[0].setSampleId(0, 0);
    families[0].setSampleId(1, 1);
    families[0].setSampleId(2, 2);
    families[1] = new Family("2", "3", "4");
    families[1].setSampleId(0, 2);
    families[1].setSampleId(1, 3);
    families[1].setSampleId(2, 4);

    families[2] = new Family("5", "6", "7");
    families[2].setSampleId(0, 5);
    families[2].setSampleId(1, 6);
    families[2].setSampleId(2, 7);
    final FamilyCaller familyCaller = new FamilyCaller(vParams, families);
    final PopulationCaller pc = new PopulationCaller(vParams, familyCaller);
    final int refNt = 2;
    final List<ModelInterface<?>> models = new ArrayList<>();
    final HypothesesSnp haploid =  new HypothesesSnp(LogApproximatePossibility.SINGLETON, params, true, refNt - 1);
    final HypothesesSnp diploid =  new HypothesesSnp(LogApproximatePossibility.SINGLETON, params, false, refNt - 1);

    final String[] values = {
        "CCCCCCCCCCCCCCCCC"
        , "CCCCCCCCCCCCCCG"
        , "CCCCCCCCCCCCCG"
        , "CCCCCCCCCCCCCCCC"
        , "CCCCCCCCCCCCGG"
        , "GGGGGGGGGGGGGGGGG"
        , "GGGGGGGGGGGGGGGGG"
        , "GGGGGGGGGGGGGGGGG"
    };
    for (String value : values) {
      final Model<Description> model = new Model<>(diploid, new StatisticsSnp(diploid.description()), new NoAlleleBalance());
      models.add(model);
      increment(model, value, 0.001, 0.001);
    }
    final byte[] ref = new byte[21];
    ref[19] = 3;
    ref[20] = 1;
    final Variant v = pc.makeCall("foo", 20, 21, ref, models, new HaploidDiploidHypotheses<>(HypothesesNone.SINGLETON, haploid, diploid));
//    System.err.println("v = " + v);
    final String[] expected = {"C:C", "C:C", "C:C", "C:C", "C:G", "G:G", "G:G", "G:G"};
    for (int i = 0; i < v.getNumberOfSamples(); i++) {
      final VariantSample sample = v.getSample(i);
      if (sample != null) {
        assertEquals(i + ": expected<" + expected[i] + "> actual<" + sample.getName() + ">",  expected[i], sample.getName());
      }
    }
    // Non ref calls so should be interesting
    assertTrue(v.toString(), v.isInteresting());

  }
  static final String REAL_WORLD_PEDIGREE = ""
                                            + "Lyon_10\tLID57249\t0\t0\t1\t1" + LS
                                            + "Lyon_10\tLID57242\t0\t0\t2\t1" + LS
                                            + "Lyon_10\tLID57250\tLID57249\tLID57242\t2\t1" + LS
                                            + "Lyon_10\tLID57243\t0\t0\t1\t1" + LS
                                            + "Lyon_10\tLID57245\tLID57249\tLID57242\t1\t1" + LS
                                            + "Lyon_10\tLID57241\t0\t0\t1\t1" + LS
                                            + "Lyon_10\tLID57246\tLID57249\tLID57242\t2\t1" + LS
                                            + "Lyon_10\tLID57247\tLID57243\tLID57250\t1\t2" + LS
                                            + "Lyon_10\tLID57244\tLID57243\tLID57250\t1\t2" + LS
                                            + "Lyon_10\tLID57248\tLID57241\tLID57246\t1\t1" + LS;

  public void testDisagreeingPPPFallbackEMRealWorld() throws IOException, PedigreeException {
    // An attempt to replicate the initial failure. See testDisagreeingPPPFallbackEM above for explanation
    final Map<String, double[]> pedigreeToModel = new HashMap<>();
    // Order of entries "A:A", "C:C", "G:G", "T:T", "A:C", "C:G", "G:T", "A:G", "C:T", "A:T"
    pedigreeToModel.put("LID57241", new double[]{-57.119, -69.947, -12.893, -69.947, -58.502, -19.113, -19.113, -7.668, -69.947, -58.502});
    pedigreeToModel.put("LID57242", new double[]{-74.691, -99.997, -25.587, -99.997, -77.456, -34.510, -34.510, -11.969, -99.997, -77.456});
    pedigreeToModel.put("LID57243", new double[]{-32.283, -70.601, -38.374, -70.601, -36.432, -41.832, -41.832, -7.662, -70.601, -36.432});
    pedigreeToModel.put("LID57244", new double[]{-125.283, -138.441, -13.354, -138.441, -126.666, -27.154, -27.154, -15.379, -138.441, -126.666});
    pedigreeToModel.put("LID57245", new double[]{-62.709, -82.363, -19.761, -82.363, -64.785, -26.661, -26.661, -9.082, -82.363, -64.785});
    pedigreeToModel.put("LID57246", new double[]{-118.455, -118.455, -0.193, -118.455, -118.455, -13.298, -13.298, -13.298, -118.455, -118.455});
    pedigreeToModel.put("LID57247", new double[]{-77.238, -90.251, -13.081, -90.251, -78.622, -21.379, -21.379, -9.749, -90.251, -78.622});
    pedigreeToModel.put("LID57248", new double[]{-48.369, -61.316, -13.212, -61.316, -49.752, -18.666, -18.666, -7.103, -61.316, -49.752});
    pedigreeToModel.put("LID57249", new double[]{-97.630, -110.742, -13.228, -110.742, -99.013, -23.589, -23.589, -11.860, -110.742, -99.013});
    pedigreeToModel.put("LID57250", new double[]{-107.690, -120.947, -13.558, -120.947, -109.074, -25.238, -25.238, -13.365, -120.947, -109.074});

    final GenomeRelationships genomeRelationships = GenomeRelationships.loadGenomeRelationships(new BufferedReader(new StringReader(REAL_WORLD_PEDIGREE)));
    final Set<Family> familySet = Family.getFamilies(genomeRelationships, false, pedigreeToModel.keySet());

    final List<Family> orderedFamilies = MultiFamilyOrdering.orderFamiliesAndSetMates(familySet);
    final Family[] families = orderedFamilies.toArray(new Family[orderedFamilies.size()]);

    final GenomePriorParams params = new GenomePriorParamsBuilder().create();
    final VariantParams vParams = new VariantParamsBuilder()
        .callLevel(VariantOutputLevel.ALL)
        .maxAmbiguity(0.3)
        .maxCoverageFilter(new StaticThreshold(50))
        .genomePriors(params)
        .maxEmIterations(6)
        .create();
    final FamilyCaller familyCaller = new FamilyCaller(vParams, families);
    final PopulationCaller pc = new PopulationCaller(vParams, familyCaller);
    final int refNt = 3;
    final List<ModelInterface<Description>> models = new ArrayList<>();
    final HypothesesSnp haploid =  new HypothesesSnp(LogApproximatePossibility.SINGLETON, params, true, refNt - 1);
    final HypothesesSnp diploid =  new HypothesesSnp(LogApproximatePossibility.SINGLETON, params, false, refNt - 1);

    for (int i = 0; i < pedigreeToModel.size(); i++) {
      models.add(null);
    }

    final List<String> calledGenomes = new ArrayList<>(pedigreeToModel.keySet());
    for (Family family : families) {
      family.setSampleIds(calledGenomes);
      final String[] members = family.getMembers();
      final int[] sampleIds = family.getSampleIds();
      for (int i = 0; i < sampleIds.length; i++) {
        final double[] posteriors = pedigreeToModel.get(members[i]);
        final double[] probSpace = new double[posteriors.length];
        for (int j = 0; j < posteriors.length; j++) {
          probSpace[j] = Math.exp(posteriors[j]);
        }
        models.set(sampleIds[i], new AbstractFamilyPosteriorTest.MockModel(diploid, probSpace));
      }
    }

    final byte[] ref = new byte[21];
    ref[19] = 3;
    ref[20] = 1;
    final Variant v = pc.makeCall("foo", 20, 21, ref, new ArrayList<ModelInterface<?>>(models), new HaploidDiploidHypotheses<>(HypothesesNone.SINGLETON, haploid, diploid));
    final String[] expected = {"A:G", "A:G", "A:G", "A:G", "G:G", "A:G", "G:G", "A:G", "A:G", "A:G"};
    for (int i = 0; i < v.getNumberOfSamples(); i++) {
      final VariantSample sample = v.getSample(i);
      if (sample != null) {
        assertEquals(i + ": expected<" + expected[i] + "> actual<" + sample.getName() + "> Variant: " + v,  expected[i], sample.getName());
      }
    }
    // Due to incomplete mock model implementation the isInteresting is incorrect
  }
}
