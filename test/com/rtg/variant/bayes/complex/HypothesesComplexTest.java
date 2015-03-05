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

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.TreeSet;

import com.rtg.mode.DNA;
import com.rtg.sam.SamUtils;
import com.rtg.util.InvalidParamsException;
import com.rtg.util.Pair;
import com.rtg.util.PortableRandom;
import com.rtg.util.StringUtils;
import com.rtg.util.Utils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.test.NanoRegression;
import com.rtg.variant.GenomePriorParams;
import com.rtg.variant.GenomePriorParamsBuilder;
import com.rtg.variant.StaticThreshold;
import com.rtg.variant.VariantAlignmentRecord;
import com.rtg.variant.VariantOutputLevel;
import com.rtg.variant.VariantParams;
import com.rtg.variant.VariantParamsBuilder;
import com.rtg.variant.bayes.multisample.population.SiteSpecificPriors;
import com.rtg.variant.match.AlignmentMatch;
import com.rtg.variant.match.Match;
import com.rtg.variant.util.VariantUtils;
import com.rtg.variant.util.arithmetic.LogPossibility;
import com.rtg.variant.util.arithmetic.SimplePossibility;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;

import junit.framework.TestCase;

/**
 */
public class HypothesesComplexTest extends TestCase {

  static final ComplexTemplate COMPLEX_TEMPLATE = new ComplexTemplate(DNA.stringDNAtoByte("CCCCCGGGGG"), "", 5, 5);

  private NanoRegression mNano = null;

  @Override
  protected void setUp() {
    Diagnostic.setLogStream();
    mNano = new NanoRegression(this.getClass(), false);
  }

  @Override
  public void tearDown() throws Exception {
    try {
      mNano.finish();
    } finally {
      mNano = null;
    }
  }

  public static VariantParams getVariantParams(final double homoSnpPrior, final double heteroSnpPrior, final double indelPrior) throws InvalidParamsException {
    final GenomePriorParams vpb = new GenomePriorParamsBuilder()
    .genomeSnpRateHomo(homoSnpPrior)
    .genomeSnpRateHetero(heteroSnpPrior)
    .genomeIndelEventRate(indelPrior / 2.0)
    .create();
    return new VariantParamsBuilder().genomePriors(vpb).defaultQuality(20).maxAmbiguity(1.0).maxCoverageFilter(new StaticThreshold(100)).create();
  }

  // Get a scorer with all the priors explicitly set
  public static HypothesesComplex getSpecifiedScorer(final ArrayList<AlignmentMatch> ml, boolean haploid, SiteSpecificPriors ssp, ComplexTemplate cot, boolean prune) {
    final VariantParams vp = getSpecifiedParams(prune);
    return HypothesesComplex.makeComplexHypotheses(cot, ml, LogPossibility.SINGLETON, haploid, vp, ssp);
  }

  // Get a scorer with all the priors explicitly set
  public static HypothesesComplex getSpecifiedScorerProb(final ArrayList<AlignmentMatch> ml, boolean haploid, SiteSpecificPriors ssp, ComplexTemplate cot, boolean prune) {
    final VariantParams vp = getSpecifiedParams(prune);
    return HypothesesComplex.makeComplexHypotheses(cot, ml, SimplePossibility.SINGLETON, haploid, vp, ssp);
  }

  private static VariantParams getSpecifiedParams(boolean prune) {
    final GenomePriorParams vpb = new GenomePriorParamsBuilder()
    .genomeIndelEventFraction(0.7)
    .genomeIndelEventRate(0.025)
    .genomeIndelDistribution(new double[] {0.5, 0.25, 0.25})
    .create();
    return new VariantParamsBuilder().genomePriors(vpb).defaultQuality(10).pruneHypotheses(prune).create();
  }

  static AlignmentMatch matchi(final int qDef) {
    final SAMRecord sam = new SAMRecord(new SAMFileHeader());
    sam.setAlignmentStart(1);
    sam.setReadString("CCCCCGGGGG");
    sam.setCigarString("10=");
    return new AlignmentMatch(new VariantAlignmentRecord(sam), "", null, qDef, 0, 0, 10);
  }

  public static AlignmentMatch match(final String ins, int readMapq) {
    final int length = ins.length();
    final SAMRecord sam = new SAMRecord(new SAMFileHeader());
    sam.setAlignmentStart(1);
    sam.setReadString("CCCCC" + ins + "GGGGG");
    sam.setCigarString("5=" + length + "I5=");
    final StringBuilder sb = new StringBuilder();
    for (int i = 0; i < length; i++) {
      sb.append('`');
    }
    return new AlignmentMatch(new VariantAlignmentRecord(sam), ins, sb.toString(), 0, 0, length, readMapq);
  }

  static class Hyp extends Pair<String, Double> {
    public Hyp(String a, Double b) {
      super(a, b);
    }
  }

  void checkHypothesis(String id, HypothesesComplex complex) throws IOException {
    StringBuilder sb = new StringBuilder();
    sb.append(complex.code().size()).append(StringUtils.LS);
    for (int i = 0; i < complex.code().size(); i++) {
      sb.append("new Hyp(\"").append(complex.name(i)).append("\", ").append(Utils.realFormat(complex.p(i), 3)).append("),").append(StringUtils.LS);
    }
    mNano.check("HypothesesComplex-" + id, sb.toString());
  }

  public void testEmptyMatches() throws Exception {
    final ArrayList<AlignmentMatch> ml = new ArrayList<>();
    final VariantParams vp = getVariantParams(0.5, 0.5, 0.1);
    final ComplexTemplate cot = new ComplexTemplate(new byte[] {}, "", 0, 0);
    final HypothesesComplex dis = HypothesesComplex.makeComplexHypotheses(cot, ml, LogPossibility.SINGLETON, true, vp, null);
    checkHypothesis("EmptyMatches", dis);
    assertEquals(0, dis.reference());
  }

  public void testAddingN() throws Exception {
    final ArrayList<AlignmentMatch> ml = new ArrayList<>();
    ml.add(match("N", 20));
    final VariantParams vp = getVariantParams(0.5, 0.5, 0.1);
    final ComplexTemplate cot = new ComplexTemplate(new byte[] {}, "", 0, 0);
    final HypothesesComplex dis = HypothesesComplex.makeComplexHypotheses(cot, ml, LogPossibility.SINGLETON, true, vp, null);
    checkHypothesis("AddingN", dis);
  }

  public void testEmptyMatchesAll() throws Exception {
    final ArrayList<AlignmentMatch> ml = new ArrayList<>();
    final GenomePriorParams vpb = new GenomePriorParamsBuilder()
    .genomeSnpRateHomo(0.5)
    .genomeSnpRateHetero(0.5)
    .genomeIndelEventRate(0.1 / 2.0)
    .create();
    final VariantParams vp = new VariantParamsBuilder().genomePriors(vpb).defaultQuality(20).callLevel(VariantOutputLevel.ALL).create();
    final ComplexTemplate cot = new ComplexTemplate(new byte[] {}, "", 0, 0);
    final HypothesesComplex dis = HypothesesComplex.makeComplexHypotheses(cot, ml, LogPossibility.SINGLETON, false, vp, null);
    checkHypothesis("EmptyMatchesAll", dis);
  }

  public void testPriorsBug() throws Exception {
    final ArrayList<AlignmentMatch> ml = new ArrayList<>();
    final GenomePriorParams gppb = new GenomePriorParamsBuilder().create();
    final VariantParamsBuilder vpb = new VariantParamsBuilder();
    vpb.mapped(new TreeSet<File>());
    final VariantParams vp = vpb.genomePriors(gppb).callLevel(VariantOutputLevel.ALL).create();
    final ComplexTemplate cot = new ComplexTemplate(new byte[] {}, "", 0, 0);
    final HypothesesComplex dis = HypothesesComplex.makeComplexHypotheses(cot, ml, LogPossibility.SINGLETON, false, vp, null);
    checkHypothesis("PriorsBug", dis);
  }


  public void testMatchesAllWithSecondOutput() throws Exception {
    final AlignmentMatch mA = match("A", 20);
    final AlignmentMatch mTT = match("TT", 20);

    final ArrayList<AlignmentMatch> ml = new ArrayList<>();
    for (int k = 0; k < 20; k++) {
      ml.add(mA);
      ml.add(mTT);
    }
    final GenomePriorParams vpb = new GenomePriorParamsBuilder().genomeSnpRateHomo(0.5).genomeSnpRateHetero(0.5)
        .genomeIndelEventRate(0.1 / 2.0)
        .create();
    final VariantParams vp = new VariantParamsBuilder().genomePriors(vpb).defaultQuality(20).callLevel(VariantOutputLevel.ALL).create();
    final HypothesesComplex dis = HypothesesComplex.makeComplexHypotheses(COMPLEX_TEMPLATE, ml, LogPossibility.SINGLETON, false, vp, null);
    checkHypothesis("MatchesAllWithSecondOutput", dis);
  }

  /*
  public void testNoInsertionMatches() {
    //since we have decided that only having one hypothesis most likely means that the complex caller was unable to create
    //alternatives due to the way reads have mapped, then this situation where all reads agree with the template causes the
    //same "unable to call" situation. In real runs the complex caller should not be invoked in this situation.
    final ArrayList<Match> ml = new ArrayList<>();
    for (int k = 0; k < 10; k++) {
      ml.add(MATCH_I);
    }
    final GenomePriorParams vpb = new GenomePriorParamsBuilder()
    .genomeSnpRateHomo(0.5)
    .genomeSnpRateHetero(0.5)
    .genomeIndelBaseRate(0.1 / 2.0)
    .genomeIndelEventRate(0.1 / 2.0)
    .create();
    final VariantParams vp = new VariantParamsBuilder().genomePriors(vpb).defaultQuality(20).callLevel(VariantOutputLevel.ALL).create();
    HypothesesComplex.makeComplexHypotheses(LogPossibility.SINGLETON, false, ml, COMPLEX_TEMPLATE, vp, new InsertionPrior(vp.priors()));
    assertEquals(1.0.totalCounts().correction(), 1E-5);
    assertEquals(10.totalCounts().count());
    assertEquals("\ti\t10\t1.000".hypothesisCorrection().output(true));
  }
   */

  public void testObviousSingleNtInsertionMatchesHomoCall() throws Exception {
    final SAMRecord samA = new SAMRecord(null);
    samA.setAlignmentStart(1);
    samA.setReadString("CCCCCAGGGGG");
    samA.setCigarString("5=1I5=");
    final AlignmentMatch mA = match("A", 20);

    final ArrayList<AlignmentMatch> ml = new ArrayList<>();
    for (int k = 0; k < 20; k++) {
      ml.add(mA);
    }
    final VariantParams vp = getVariantParams(0.5, 0.5, 0.1);
    final HypothesesComplex dis = HypothesesComplex.makeComplexHypotheses(COMPLEX_TEMPLATE, ml, LogPossibility.SINGLETON, false, vp, null);

    checkHypothesis("ObviousSingleNtInsertionMatchesHomoCall", dis);
  }

  public void testObviousSingleNtInsertionMatchesHeteroCall() throws Exception {
    final AlignmentMatch mA = match("A", 20);
    final AlignmentMatch mTT = match("TT", 20);
    final ArrayList<AlignmentMatch> ml = new ArrayList<>();
    for (int k = 0; k < 20; k++) {
      ml.add(mA);
      ml.add(mTT);
    }
    final VariantParams vp = getVariantParams(0.5, 0.5, 0.1);

    final ComplexTemplate cot = new ComplexTemplate(DNA.stringDNAtoByte("GGATCGGGGG"), "", 5, 5);
    final HypothesesComplex dis = HypothesesComplex.makeComplexHypotheses(cot, ml, LogPossibility.SINGLETON, false, vp, null);

    checkHypothesis("ObviousSingleNtInsertionMatchesHeteroCall", dis);
  }

  public void testObviousSingleNtReferenceId() throws Exception {
    final AlignmentMatch mA = match("A", 20);
    final ArrayList<AlignmentMatch> ml = new ArrayList<>();
    ml.add(mA);
    final VariantParams vp = getVariantParams(0.5, 0.5, 0.1);
    final ComplexTemplate cot = new ComplexTemplate(DNA.stringDNAtoByte("GGATCGGGGG"), "", 5, 6);
    final HypothesesComplex dis = HypothesesComplex.makeComplexHypotheses(cot, ml, LogPossibility.SINGLETON, false, vp, null);
    assertEquals(1, dis.reference());
  }

  private static final PortableRandom RANDOM = new PortableRandom();

  private AlignmentMatch randomMatch() {
    if (RANDOM.nextInt(10) == 0) {
      //return InsertionMatch.defaultInsertion(RANDOM.nextInt(64));
      final SAMRecord sam = new SAMRecord(new SAMFileHeader());
      sam.setAlignmentStart(1);
      sam.setReadString("CCCCCGGGGG");
      sam.setCigarString("10=");
      return new AlignmentMatch(new VariantAlignmentRecord(sam), "", null, RANDOM.nextInt(64), 0, 0, 10);
    }
    final StringBuilder sb = new StringBuilder();
    do {
      sb.append("ACGT".charAt(RANDOM.nextInt(4)));
    } while (RANDOM.nextBoolean());
    return match(sb.toString(), 20);
  }

  public void testRandom() throws Exception {
    for (int j = 0; j < 10; j++) {
      final VariantParams vp = getVariantParams(RANDOM.nextDouble() / 10.0, RANDOM.nextDouble() / 10.0, RANDOM.nextDouble() / 10.0);
      final ArrayList<AlignmentMatch> ml = new ArrayList<>();
      for (int k = 0; k < RANDOM.nextInt(500); k++) {
        final AlignmentMatch rm = randomMatch();
        ml.add(rm);
      }

      final HypothesesComplex dis = HypothesesComplex.makeComplexHypotheses(COMPLEX_TEMPLATE, ml, LogPossibility.SINGLETON, false, vp, null);

      for (int i = 0 ; i < dis.code().size(); i++) {
        assertFalse(Double.isInfinite(dis.p(i)));
        assertFalse(Double.isNaN(dis.p(i)));
      }
    }
  }


  // Test probability = 0.0 case
  public void testFromConvergence() throws InvalidParamsException, IOException {
    final ArrayList<AlignmentMatch> ml = new ArrayList<>();

    for (int i = 0; i < 8; i++) {
      ml.add(matchi(20));
    }
    final AlignmentMatch cc = match("AA", 20);
    for (int i = 0; i < 15; i++) {
      ml.add(cc);
    }
    ml.add(match("TGCTAA", 20));
    //System.err.println(ml);
    final GenomePriorParams vpb = new GenomePriorParamsBuilder().create();
    final VariantParams vp = new VariantParamsBuilder().genomePriors(vpb).defaultQuality(20).maxAmbiguity(100.0).maxCoverageFilter(new StaticThreshold(100)).create();
    final HypothesesComplex dis = HypothesesComplex.makeComplexHypotheses(COMPLEX_TEMPLATE, ml, LogPossibility.SINGLETON, false, vp, null);
    checkHypothesis("FromConvergence", dis);
  }

  // Test probability = 0.0 case including non-identity posterior
  public void testFromConvergence2() throws InvalidParamsException, IOException {
    final ArrayList<AlignmentMatch> ml = new ArrayList<>();

    for (int i = 0; i < 8; i++) {
      ml.add(matchi(20));
    }
    final AlignmentMatch cc = match("AA", 20);
    for (int i = 0; i < 15; i++) {
      ml.add(cc);
    }
    ml.add(match("TGCTAA", 20));

    final GenomePriorParams vpb = new GenomePriorParamsBuilder().create();
    final VariantParams vp = new VariantParamsBuilder().genomePriors(vpb).defaultQuality(20).create();
    final HypothesesComplex dis = HypothesesComplex.makeComplexHypotheses(COMPLEX_TEMPLATE, ml, LogPossibility.SINGLETON, false, vp, null);
    checkHypothesis("FromConvergence2", dis);
  }

  public void testNonDefaultEventRates() throws Exception {
    final ArrayList<AlignmentMatch> ml = new ArrayList<>();
    final HypothesesComplex dis = getSpecifiedScorer(ml, false, null, COMPLEX_TEMPLATE, false);
    checkHypothesis("NonDefaultEventRates", dis);
  }

  public void testCloneConstructor() throws Exception {
    final ArrayList<AlignmentMatch> ml = new ArrayList<>();
    final HypothesesComplex dis = getSpecifiedScorer(ml, false, null, COMPLEX_TEMPLATE, false);
    checkHypothesis("CloneConstructor", dis);

    final HypothesesComplex hc = new HypothesesComplex(dis, new double[]{0.5});
    checkHypothesis("CloneConstructor2", hc);
  }


  public void testAmbiguousRead() throws Exception {
    final SAMRecord sam = new SAMRecord(new SAMFileHeader());
    sam.setAlignmentStart(1);
    sam.setReadString("tatgaaattctgggttgaaaatttttaagaatg");
    sam.setCigarString("21=2I4N10=");
    sam.setAttribute(SamUtils.CG_SUPER_CIGAR, "5=2B18=2I4N10=");
    sam.setAttribute(SamUtils.CG_READ_DELTA, "AA");
    sam.setAttribute(SamUtils.CG_OVERLAP_QUALITY, "55");
    sam.setFlags(73);

    final AlignmentMatch ma = new AlignmentMatch(new VariantAlignmentRecord(sam), "A", "5", 0, 0, 1, 0) {
      @Override
      public double mapError() {
        return 1.0;
      }
    };

    final ArrayList<AlignmentMatch> ml = new ArrayList<>();
    ml.add(ma);
    final VariantParams vp = new VariantParamsBuilder().genomePriors(new GenomePriorParamsBuilder().create()).maxCoverageFilter(new StaticThreshold(10)).create();
    final ComplexTemplate cot = new ComplexTemplate(DNA.stringDNAtoByte("gatatgaaattctgggttgaaaattcttttaagaatgttgaatattggcccccacggtcttctggcttgtagggtttctgcagag"), "", 5, 6);
    final HypothesesComplex dis = HypothesesComplex.makeComplexHypotheses(cot, ml, LogPossibility.SINGLETON, false, vp, null);
    checkHypothesis("AmbiguousRead", dis);
  }

  public void testArithmetic() {
    final ArrayList<AlignmentMatch> ml = new ArrayList<>();
    ml.add(match("A", 20));
    ml.add(match("TT", 20));
    final GenomePriorParams vpb = new GenomePriorParamsBuilder().genomeSnpRateHomo(0.5).genomeSnpRateHetero(0.5).genomeIndelEventRate(0.1 / 2.0).create();
    final VariantParams vp = new VariantParamsBuilder().genomePriors(vpb).defaultQuality(20).callLevel(VariantOutputLevel.ALL).create();
    final HypothesesComplex dis = HypothesesComplex.makeComplexHypotheses(COMPLEX_TEMPLATE, ml, LogPossibility.SINGLETON, false, vp, null);
    final HypothesesComplex dis2 = HypothesesComplex.makeComplexHypotheses(COMPLEX_TEMPLATE, ml, SimplePossibility.SINGLETON, false, vp, null);
    for (int i = 0; i < dis.code().size(); i++) {
      assertEquals(Math.exp(dis.p(i)), dis2.p(i), 1e-11);
    }
  }

  //test reference match is always present in description
  public void testRegression() {
    final ArrayList<AlignmentMatch> ml = new ArrayList<>();
    ml.add(match("A", 40));
    final AlignmentMatch match = match("TT", 30);
    ml.add(match);
    final ComplexTemplate cot = new ComplexTemplate(DNA.stringDNAtoByte("gata"), "", 1, 2);
    final GenomePriorParams vpb = new GenomePriorParamsBuilder().genomeSnpRateHomo(0.5).genomeSnpRateHetero(0.5).genomeIndelEventRate(0.1 / 2.0).create();
    final VariantParams vp = new VariantParamsBuilder().genomePriors(vpb).defaultQuality(20).callLevel(VariantOutputLevel.ALL).create();
    final HypothesesComplex dis = HypothesesComplex.makeComplexHypotheses(cot, ml, SimplePossibility.SINGLETON, false, vp, null);
    assertEquals(2, dis.description().size());
    assertEquals(0.0, dis.description().match(0).mapError());
    assertEquals(VariantUtils.phredToProb(30), dis.description().match(1).mapError());

  }

  public void testHypothesisCut() {

    final String ref = "CCCACTGTTGC";

    final ArrayList<AlignmentMatch> ml = createMatches(ref);

    final ComplexTemplate cot = new ComplexTemplate(DNA.stringDNAtoByte(ref), "", 0, ref.length());
    final DescriptionComplex dc = HypothesesComplex.createDescription(ml, cot, true, null, Integer.MAX_VALUE);
    assertEquals(6, dc.size());
  }

  private ArrayList<AlignmentMatch> createMatches(final String ref) {
    final ArrayList<AlignmentMatch> ml = new ArrayList<>();
    for (int i = 0; i <= 934; i++) {
      ml.add(match(ref, 40));
    }
    for (int i = 0; i <= 6; i++) {
      ml.add(match("CCACTGTTGC", 40));
    }
    for (int i = 0; i <= 5; i++) {
      ml.add(match("CCCACTGCTGC", 40));
    }
    for (int i = 0; i <= 5; i++) {
      ml.add(match("CCCACTGTCGC", 40));
    }
    for (int i = 0; i <= 5; i++) {
      ml.add(match("CCCCACTGTTGC", 40));
    }
    ml.add(match("CCCACTAGTTGC", 40));
    ml.add(match("CCCACTAGTTGC", 40));

    ml.add(match("CCAACTTGTTGC", 40));
    ml.add(match("CCACATGTTGC", 40));
    ml.add(match("CCACTGTGC", 40));
    ml.add(match("CCATGTTGC", 40));
    ml.add(match("CCCACATGTTGC", 40));
    ml.add(match("CCCACCGTTGC", 40));
    ml.add(match("CCCACTATTGC", 40));
    ml.add(match("CCCACTGTGC", 40));
    ml.add(match("CCCACTGTTAC", 40));
    ml.add(match("CCCACTGTTAGC", 40));
    ml.add(match("CCCACTGTTGCA", 40));
    ml.add(match("CCCACTGTTGT", 40));
    ml.add(match("CCCACTTGTTGC", 40));
    ml.add(match("CCCACTTGTTGCA", 40));
    ml.add(match("CCCCCTGTTGC", 40));
    ml.add(match("CCCGCTGTTGC", 40));
    ml.add(match("CCCTACTGTTGC", 40));
    ml.add(match("CCTACTGTTGC", 40));
    ml.add(match("TCCACTGTTGC", 40));
    ml.add(match("TCCCACTGTTG", 40));
    return ml;
  }

  public void testHypothesisCut2() {
    final String ref = "CAAATCCCAAG";
    final ArrayList<AlignmentMatch> ml = new ArrayList<>();
    hypothesisCut2ForLoops1(ref, ml);
    hypothesisCut2ForLoops2(ml);
    ml.add(match("ACAAATCCCAA", 40));
    ml.add(match("ACAATCCCAAG", 40));
    ml.add(match("CAAACCCAAA", 40));
    ml.add(match("CAAATCCAA", 40));
    ml.add(match("CAAATCCAATC", 40));
    ml.add(match("CAAATCCCAA", 40));
    ml.add(match("CAAATCCCAGA", 40));
    ml.add(match("CAAATCCCGAA", 40));
    ml.add(match("CAAATCCTAAG", 40));
    ml.add(match("CAATCCAAA", 40));
    ml.add(match("CGAATCCTAAG", 40));
    ml.add(match("GCAAATCCCAAG", 40));
    ml.add(match("TAAATCCCAAA", 40));
    ml.add(match("TCAAATCCCAAA", 40));
    ml.add(match("TCAAATCCCAAG", 40));

    final ComplexTemplate cot = new ComplexTemplate(DNA.stringDNAtoByte(ref), "", 0, ref.length());
    final DescriptionComplex dc = HypothesesComplex.createDescription(ml, cot, true, null, 6);
    assertEquals(6, dc.size());
  }

  private void hypothesisCut2ForLoops2(ArrayList<AlignmentMatch> ml) {
    for (int i = 0; i <= 2; i++) {
      ml.add(match("CAAACCCCAAA", 40));
    }
    for (int i = 0; i <= 2; i++) {
      ml.add(match("CAAACCCCAAG", 40));
    }
    for (int i = 0; i <= 2; i++) {
      ml.add(match("CAAATCCAAA", 40));
    }
    for (int i = 0; i <= 2; i++) {
      ml.add(match("CAAATCCCAAAA", 40));
    }
    for (int i = 0; i <= 2; i++) {
      ml.add(match("CAAATCCCAAAG", 40));
    }
    for (int i = 0; i <= 2; i++) {
      ml.add(match("CAAATCCCAATG", 40));
    }
    for (int i = 0; i <= 2; i++) {
      ml.add(match("CAAATTCCAAA", 40));
    }
    for (int i = 0; i <= 2; i++) {
      ml.add(match("CAAATTCCCAAA", 40));
    }
    for (int i = 0; i <= 2; i++) {
      ml.add(match("CAAGTCCCAAG", 40));
    }
  }

  private void hypothesisCut2ForLoops1(String ref, ArrayList<AlignmentMatch> ml) {
    for (int i = 0; i <= 592; i++) {
      ml.add(match(ref, 40));
    }
    for (int i = 0; i <= 565; i++) {
      ml.add(match("CAAATCCCAAA", 40));
    }
    for (int i = 0; i <= 13; i++) {
      ml.add(match("CAATCCCAAG", 40));
    }
    for (int i = 0; i <= 8; i++) {
      ml.add(match("CAATCCCAAA", 40));
    }
    for (int i = 0; i <= 7; i++) {
      ml.add(match("CAAATCCCCAAA", 40));
    }
    for (int i = 0; i <= 4; i++) {
      ml.add(match("CAAATCCAAG", 40));
    }
    for (int i = 0; i <= 3; i++) {
      ml.add(match("CAAATCCCGAG", 40));
    }
    for (int i = 0; i <= 3; i++) {
      ml.add(match("CAAAATCCCAAA", 40));
    }
    for (int i = 0; i <= 2; i++) {
      ml.add(match("ACAAATCCCAAA", 40));
    }
    for (int i = 0; i <= 2; i++) {
      ml.add(match("CAAAATCCCAAG", 40));
    }
  }

  public void testIndelAndUnknownRefBug1591() {
    final ComplexTemplate cot = new ComplexTemplate(new byte[1000], "ref", 468, 510);
    final Match m = new AlignmentMatch(null, "AAAAA", null, 20, 0, 5, 0);
    final DescriptionComplex desc = new DescriptionComplex(Collections.singletonList(m));
    assertNotNull(HypothesesComplex.makePriorsAllPaths(desc, true, cot, SimplePossibility.SINGLETON, -1, new GenomePriorParamsBuilder().create()));
  }

  // TODO Following test should be enabled after fixing TODO in HypothesesComplex.createDescription method
//  public void testHypothesisCutWithSsp() {
//    final String ref = "CCCACTGTTGC";
//
//    final ArrayList<AlignmentMatch> ml = createMatches(ref);
//
//    final ComplexTemplate cot = new ComplexTemplate(DNA.stringDNAtoByte(ref), "", 0, ref.length());
//
//    final SiteSpecificPriors ssp = new SiteSpecificPriors() {
//      @Override
//      public HaploidDiploidHypotheses<?> getHypotheses(String templateName, int zeroPos) {
//        return null;
//      }
//
//      @Override
//      public List<AlleleCounts> getCounts(String templateName, int start, int end) {
//        final HashMap<String, Integer> counts = new HashMap<>();
//        counts.put("CCAACTTGTTGC", 12);
//        counts.put("CCACATGTTGC", 11);
//        counts.put("AAAAAAAAAA", 10);
//        counts.put("T", 1);
//        final ArrayList<AlleleCounts> list = new ArrayList<>();
//        list.add(new AlleleCounts(1, counts, 5, false));
//        return list;
//      }
//
//    };
//
//    final DescriptionComplex dc = HypothesesComplex.createDescription(ml, cot, true, ssp);
//
//    assertEquals(10, dc.size());
//    final HashSet<String> names = new HashSet<>();
//    names.add("CAAAAAAAAAAGTTGC");
//    names.add("CCACTGTTGC");
//    names.add("CCCAACTTGTTGCGTTGC");
//    names.add("CCCACATGTTGCGTTGC");
//    names.add("CCCACTAGTTGC");
//    names.add("CCCACTGCTGC");
//    names.add("CCCACTGTCGC");
//    names.add("CCCACTGTTGC");
//    names.add("CCCCACTGTTGC");
//    names.add("CTGTTGC");
//
//    for (int i = 0; i < dc.size(); i++) {
//      assertTrue(names.contains(dc.name(i)));
//    }
//  }
//
//  public void testHypothesisCutWithSspNoPrune() {
//    final String ref = "CCCACTGTTGC";
//
//    final ArrayList<AlignmentMatch> ml = new ArrayList<>();
//    for (int i = 0; i < 100; i++) {
//      ml.add(match(ref, 40));
//    }
//    ml.add(match("CCACTGTTGC", 40));
//    ml.add(match("CCCACTGCTGC", 40));
//    ml.add(match("CCCACTGTCGC", 40));
//    ml.add(match("CCCCACTGTTGC", 40));
//    ml.add(match("CCCACTAGTTGC", 40));
//    ml.add(match("CCCACTAGTTGC", 40));
//
//    final ComplexTemplate cot = new ComplexTemplate(DNA.stringDNAtoByte(ref), "", 0, ref.length());
//
//    final SiteSpecificPriors ssp = new SiteSpecificPriors() {
//      @Override
//      public HaploidDiploidHypotheses<?> getHypotheses(String templateName, int zeroPos) {
//        return null;
//      }
//
//      @Override
//      public List<AlleleCounts> getCounts(String templateName, int start, int end) {
//        final HashMap<String, Integer> counts = new HashMap<>();
//        counts.put("CCAACTTGTTGC", 12);
//        counts.put("CCACATGTTGC", 11);
//        counts.put("AAAAAAAAAA", 10);
//        counts.put("T", 1);
//        final ArrayList<AlleleCounts> list = new ArrayList<>();
//        list.add(new AlleleCounts(1, counts, 5, false));
//        return list;
//      }
//
//    };
//
//    final DescriptionComplex dc = HypothesesComplex.createDescription(ml, cot, false, ssp);
//
//    assertEquals(10, dc.size());
//    final HashSet<String> names = new HashSet<>();
//    names.add("CAAAAAAAAAAGTTGC");
//    names.add("CCACTGTTGC");
//    names.add("CCCAACTTGTTGCGTTGC");
//    names.add("CCCACATGTTGCGTTGC");
//    names.add("CCCACTAGTTGC");
//    names.add("CCCACTGCTGC");
//    names.add("CCCACTGTCGC");
//    names.add("CCCACTGTTGC");
//    names.add("CCCCACTGTTGC");
//    names.add("CTGTTGC");
//
//    for (int i = 0; i < dc.size(); i++) {
//      assertTrue(names.contains(dc.name(i)));
//    }
//  }
//

//  public void testRefSeqFromAllele() {
//    final String ref = "CAAATCCCAAG";
//    final ComplexTemplate cot = new ComplexTemplate(DNA.stringDNAtoByte(ref), "", 0, ref.length());
//
//    byte[] res = HypothesesComplex.refSeqFromAllele(cot, 1, 1, "T");
//    assertNotNull(res);
//    assertEquals(ref.length(), res.length);
//    for (int i = 0; i < ref.length(); i++) {
//      if (i == 1) {
//        assertEquals(DNA.getDNA('T'), res[i]);
//      } else {
//        assertEquals(DNA.getDNA(ref.charAt(i)), res[i]);
//      }
//    }
//
//    res = HypothesesComplex.refSeqFromAllele(cot, 1, 1, "A");
//    assertNotNull(res);
//    assertEquals(ref.length(), res.length);
//    for (int i = 0; i < ref.length(); i++) {
//      assertEquals(DNA.getDNA(ref.charAt(i)), res[i]);
//    }
//
//    res = HypothesesComplex.refSeqFromAllele(cot, 1, 1, "AC");
//    assertNotNull(res);
//    assertEquals(ref.length() + 1, res.length);
//    for (int i = 0; i < res.length; i++) {
//      if (i == 2) {
//        assertEquals(DNA.getDNA('C'), res[i]);
//      } else if (i < 2) {
//        assertEquals(DNA.getDNA(ref.charAt(i)), res[i]);
//      } else {
//        assertEquals("" + i + " " + DnaUtils.bytesToSequenceIncCG(res), DNA.getDNA(ref.charAt(i - 1)), res[i]);
//      }
//    }
//
//    res = HypothesesComplex.refSeqFromAllele(cot, 1, 1, "");
//    assertNotNull(res);
//    assertEquals(ref.length() - 1, res.length);
//    for (int i = 0; i < res.length; i++) {
//      if (i < 1) {
//        assertEquals(DNA.getDNA(ref.charAt(i)), res[i]);
//      } else {
//        assertEquals("" + i + " " + DnaUtils.bytesToSequenceIncCG(res), DNA.getDNA(ref.charAt(i + 1)), res[i]);
//      }
//    }
//  }
//  public void testCreateDescriptionWithSsp() {
//    final String ref = "CAAATCCCAAG";
//    final ComplexTemplate cot = new ComplexTemplate(DNA.stringDNAtoByte(ref), "", 0, ref.length());
//
//    final SiteSpecificPriors ssp = new SiteSpecificPriors() {
//      @Override
//      public HaploidDiploidHypotheses<?> getHypotheses(String templateName, int zeroPos) {
//        return null;
//      }
//
//      @Override
//      public List<AlleleCounts> getCounts(String templateName, int start, int end) {
//        final HashMap<String, Integer> counts = new HashMap<>();
//        counts.put("T", 1);
//        counts.put("A", 1);
//        counts.put("AC", 1);
//        counts.put("", 1);
//        final ArrayList<AlleleCounts> list = new ArrayList<>();
//        list.add(new AlleleCounts(1, counts, 1, false));
//        return list;
//      }
//    };
//
//    final DescriptionComplex description = HypothesesComplex.createDescription(new ArrayList<AlignmentMatch>(), cot, false, ssp);
//    assertNotNull(description);
//    assertEquals(4, description.size());
//    final HashSet<String> names = new HashSet<>();
//    for (int i = 0; i < description.size(); i++) {
//      names.add(description.name(i));
//    }
//    assertEquals(description.size(), names.size());
//
//    assertTrue(names.contains("CAAATCCCAAG"));
//    assertTrue(names.contains("CTAATCCCAAG"));
//    assertTrue(names.contains("CACAATCCCAAG"));
//    assertTrue(names.contains("CAATCCCAAG"));
//  }
}
