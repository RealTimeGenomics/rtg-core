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
import java.util.List;
import java.util.TreeSet;

import com.rtg.mode.DNA;
import com.rtg.mode.DnaUtils;
import com.rtg.sam.SamUtils;
import com.rtg.util.InvalidParamsException;
import com.rtg.util.PortableRandom;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.variant.DefaultMachineErrorChooser;
import com.rtg.variant.GenomePriorParams;
import com.rtg.variant.GenomePriorParamsBuilder;
import com.rtg.variant.MachineErrorChooserInterface;
import com.rtg.variant.MachineErrorParams;
import com.rtg.variant.StaticThreshold;
import com.rtg.variant.VariantAlignmentRecord;
import com.rtg.variant.VariantOutputLevel;
import com.rtg.variant.VariantParams;
import com.rtg.variant.VariantParamsBuilder;
import com.rtg.variant.bayes.EvidenceInterface;
import com.rtg.variant.bayes.MockHypotheses;
import com.rtg.variant.bayes.Model;
import com.rtg.variant.bayes.ModelInterface;
import com.rtg.variant.bayes.NoAlleleBalance;
import com.rtg.variant.bayes.complex.HypothesesComplexTest.Hyp;
import com.rtg.variant.match.AlignmentMatch;
import com.rtg.variant.match.AlleleAsReadMatch;
import com.rtg.variant.match.Match;
import com.rtg.variant.util.arithmetic.LogPossibility;
import com.rtg.variant.util.arithmetic.PossibilityArithmetic;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import junit.framework.TestCase;

/**
 */
public class EvidenceComplexTest extends TestCase {

  @Override
  public void setUp() {
    Diagnostic.setLogStream();
  }

  private static MachineErrorParams getErrors() {
    final MachineErrorParams me =  MachineErrorParams.builder().errorDelEventRate(0.009).errorInsEventRate(0.006).create();
    assertEquals(0.015, me.errorDelBaseRate(), 1e-5);
    assertEquals(0.01, me.errorInsBaseRate(), 1e-5);
    return me;
  }

  private static MachineErrorParams getMachineErrors() {
    return getErrors();
  }

  /**
   * @return test machine error chooser
   */
  public static MachineErrorChooserInterface getChooser() {
    return new DefaultMachineErrorChooser(getMachineErrors());
  }

  private ModelInterface<?> getModel(final HypothesesComplex hyp, final StatisticsComplex statistics) {
    return new Model<>(hyp, statistics, new NoAlleleBalance());
  }

  private static void checkPosterior(ModelInterface<?> m, Hyp... hyp) {
    for (int i = 0 ; i < m.size(); i++) {
      final String msg = i + ":" + m.name(i);
      //System.err.println(msg + ":" + m.posteriorLn(i));
      assertEquals(msg, hyp[i].getA(), m.name(i));
      //assertEquals(msg, hyp[i].getB(), m.posteriorLn(i), 1e-3);
    }
  }

  public void testMatchesAllWithSecondOutput() {
    final AlignmentMatch mA = HypothesesComplexTest.match("A", 5);
    final AlignmentMatch mTT = HypothesesComplexTest.match("TT", 5);

    final ArrayList<AlignmentMatch> ml = new ArrayList<>();
    for (int i = 0; i < 20; i++) {
      ml.add(mA);
      ml.add(mTT);
    }
    final MachineErrorParams me = getErrors();
    final MachineErrorChooserInterface ch = new DefaultMachineErrorChooser(me);
    final GenomePriorParams vpb = new GenomePriorParamsBuilder().genomeSnpRateHomo(0.5).genomeSnpRateHetero(0.5)
        .genomeIndelEventRate(0.1 / 2.0)
        .create();
    final VariantParams vp = new VariantParamsBuilder().genomePriors(vpb).defaultQuality(20).callLevel(VariantOutputLevel.ALL).create();
    final HypothesesComplex hyp = HypothesesComplex.makeComplexHypotheses(HypothesesComplexTest.COMPLEX_TEMPLATE, ml, LogPossibility.SINGLETON, false, vp, null);
    final StatisticsComplex statistics = new StatisticsComplex(hyp.description(), HypothesesComplexTest.COMPLEX_TEMPLATE.getLength());
    final ModelInterface<?> m = getModel(hyp, statistics);
    double sum = 0.0;
    for (final AlignmentMatch match : ml) {
      final EvidenceComplex evidence = new EvidenceComplex(hyp, match, HypothesesComplexTest.COMPLEX_TEMPLATE, vp, ch);
      //System.err.println(evidence.sumLn());
      sum += evidence.sumLn();
      m.increment(evidence);
    }
    assertEquals(-16.1813, sum, 1e-3);
    checkPosterior(m, new Hyp(":", -259.067),
        new Hyp("A:A", -180.472),
        new Hyp("TT:TT", -129.989),
        new Hyp(":A", -191.105),
        new Hyp("A:TT", -54.033),
        new Hyp(":TT", -135.688)
        );
  }

  //with new P(E) calculation
  public void testMatchesAllWithSecondOutputPE() {
    System.getProperties().setProperty("complex_pe", "true");
    final AlignmentMatch mA = HypothesesComplexTest.match("A", 5);
    final AlignmentMatch mTT = HypothesesComplexTest.match("TT", 5);

    final ArrayList<AlignmentMatch> ml = new ArrayList<>();
    for (int i = 0; i < 20; i++) {
      ml.add(mA);
      ml.add(mTT);
    }
    final MachineErrorParams machineErrors = getErrors();
    final MachineErrorChooserInterface ch = new DefaultMachineErrorChooser(machineErrors);
    final GenomePriorParams vpb = new GenomePriorParamsBuilder().genomeSnpRateHomo(0.5).genomeSnpRateHetero(0.5)
        .genomeIndelEventRate(0.1 / 2.0)
        .create();
    final VariantParams vp = new VariantParamsBuilder().genomePriors(vpb).defaultQuality(20).callLevel(VariantOutputLevel.ALL).create();
    final HypothesesComplex hypComplex = HypothesesComplex.makeComplexHypotheses(HypothesesComplexTest.COMPLEX_TEMPLATE, ml, LogPossibility.SINGLETON, false, vp, null);
    final StatisticsComplex stat = new StatisticsComplex(hypComplex.description(), HypothesesComplexTest.COMPLEX_TEMPLATE.getLength());
    final ModelInterface<?> m = getModel(hypComplex, stat);
    double sum = 0.0;
    for (final AlignmentMatch match : ml) {
      final EvidenceComplex evidence = new EvidenceComplex(hypComplex, match, HypothesesComplexTest.COMPLEX_TEMPLATE, vp, ch);
      sum += evidence.sumLn();
      m.increment(evidence);
    }
    assertEquals(-16.18125, sum, 1e-3);
    checkPosterior(m, new Hyp(":", -259.067),
        new Hyp("A:A", -180.472),
        new Hyp("TT:TT", -129.989),
        new Hyp(":A", -191.105),
        new Hyp("A:TT", -54.033),
        new Hyp(":TT", -135.688)
        );
    System.getProperties().setProperty("complex_pe", "false");

  }

  public void testObviousSingleNtInsertionMatchesHomoCall() throws Exception {
    final SAMRecord samA = new SAMRecord(null);
    samA.setAlignmentStart(1);
    samA.setReadString("CCCCCAGGGGG");
    samA.setCigarString("5=1I5=");
    final AlignmentMatch mA = HypothesesComplexTest.match("A", 20);

    final ArrayList<AlignmentMatch> ml = new ArrayList<>();
    for (int k = 0; k < 20; k++) {
      ml.add(mA);
    }
    final VariantParams vp = HypothesesComplexTest.getVariantParams(0.5, 0.5, 0.1);
    final HypothesesComplex hyp = HypothesesComplex.makeComplexHypotheses(HypothesesComplexTest.COMPLEX_TEMPLATE, ml, LogPossibility.SINGLETON, false, vp, null);
    final StatisticsComplex statistics = new StatisticsComplex(hyp.description(), HypothesesComplexTest.COMPLEX_TEMPLATE.getLength());
    final ModelInterface<?> m = getModel(hyp, statistics);
    double sum = 0.0;
    for (final AlignmentMatch match : ml) {
      final EvidenceComplex evidence = new EvidenceComplex(hyp, match, HypothesesComplexTest.COMPLEX_TEMPLATE, vp, getChooser());
      sum += evidence.sumLn();
      m.increment(evidence);
    }
    assertEquals(-7.66098, sum, 0.0001);
    checkPosterior(m, new Hyp(":", -103.0518),
        new Hyp("A:A", -5.735),
        new Hyp(":A", -20.544));
  }

  public void testEmptyMatches() throws Exception {
    final ArrayList<AlignmentMatch> ml = new ArrayList<>();
    final VariantParams vp = HypothesesComplexTest.getVariantParams(0.5, 0.5, 0.1);
    final ComplexTemplate cot = new ComplexTemplate(new byte[] {}, "", 0, 0);
    final HypothesesComplex hyp = HypothesesComplex.makeComplexHypotheses(cot, ml, LogPossibility.SINGLETON, true, vp, null);
    final StatisticsComplex statistics = new StatisticsComplex(hyp.description(), cot.getLength());
    double sum = 0.0;
    final ModelInterface<?> m = getModel(hyp, statistics);
    for (final AlignmentMatch match : ml) {
      final EvidenceComplex evidence = new EvidenceComplex(hyp, match, cot, vp, getChooser());
      sum += evidence.sumLn();
      m.increment(evidence);
    }
    checkPosterior(m, new Hyp("", -0.051 - sum));
  }

  public void testAddingN() throws Exception {
    final ArrayList<AlignmentMatch> ml = new ArrayList<>();
    ml.add(HypothesesComplexTest.match("N", 20));
    final VariantParams vp = HypothesesComplexTest.getVariantParams(0.5, 0.5, 0.1);
    final ComplexTemplate cot = new ComplexTemplate(new byte[] {}, "", 0, 0);
    final HypothesesComplex hyp = HypothesesComplex.makeComplexHypotheses(cot, ml, LogPossibility.SINGLETON, true, vp, null);
    final StatisticsComplex statistics = new StatisticsComplex(hyp.description(), cot.getLength());
    final ModelInterface<?> m = getModel(hyp, statistics);
    double sum = 0.0;
    for (final AlignmentMatch match : ml) {
      final EvidenceComplex evidence = new EvidenceComplex(hyp, match, cot, vp, getChooser());
      sum += evidence.sumLn();
      m.increment(evidence);
    }
    assertEquals(-12.564, sum, 0.001);
    checkPosterior(m, new Hyp("", -0.0518));
  }

  public void testEmptyMatchesAll() {
    final ArrayList<AlignmentMatch> ml = new ArrayList<>();
    final GenomePriorParams vpb = new GenomePriorParamsBuilder()
    .genomeSnpRateHomo(0.5)
    .genomeSnpRateHetero(0.5)
    .genomeIndelEventRate(0.1 / 2.0)
    .create();
    final MachineErrorParams me = getMachineErrors();
    final MachineErrorChooserInterface ch = new DefaultMachineErrorChooser(me);
    final VariantParams vp = new VariantParamsBuilder().genomePriors(vpb).defaultQuality(20).callLevel(VariantOutputLevel.ALL).create();
    final ComplexTemplate cot = new ComplexTemplate(new byte[] {}, "", 0, 0);
    final HypothesesComplex hyp = HypothesesComplex.makeComplexHypotheses(cot, ml, LogPossibility.SINGLETON, false, vp, null);
    final StatisticsComplex statistics = new StatisticsComplex(hyp.description(), cot.getLength());
    final ModelInterface<?> m = getModel(hyp, statistics);
    for (final AlignmentMatch match : ml) {
      m.increment(new EvidenceComplex(hyp, match, cot, vp, ch));
    }
    checkPosterior(m, new Hyp(":", -0.051));
  }

  public void testPriorsBug() {
    final ArrayList<AlignmentMatch> ml = new ArrayList<>();
    final GenomePriorParams gppb = new GenomePriorParamsBuilder().create();
    final VariantParamsBuilder vpb = new VariantParamsBuilder();
    vpb.mapped(new TreeSet<File>());
    final VariantParams vp = vpb.genomePriors(gppb).callLevel(VariantOutputLevel.ALL).create();
    final ComplexTemplate cot = new ComplexTemplate(new byte[] {}, "", 0, 0);
    final HypothesesComplex hyp = HypothesesComplex.makeComplexHypotheses(cot, ml, LogPossibility.SINGLETON, true, vp, null);
    final StatisticsComplex statistics = new StatisticsComplex(hyp.description(), cot.getLength());
    final ModelInterface<?> m = getModel(hyp, statistics);
    for (final AlignmentMatch match : ml) {
      m.increment(new EvidenceComplex(hyp, match, HypothesesComplexTest.COMPLEX_TEMPLATE, vp, getChooser()));
    }
    checkPosterior(m, new Hyp("", 0.000));

  }

  public void testReadHypothesis() throws Exception {

    //keeps in treeset
    // "" - ref
    // A
    // CC
    // GG
    // TT
    final AlignmentMatch mA = HypothesesComplexTest.match("A", 20);
    final AlignmentMatch mTT = HypothesesComplexTest.match("TT", 20);
    final AlignmentMatch mGG = HypothesesComplexTest.match("GG", 20);
    final AlignmentMatch mCC = HypothesesComplexTest.match("CC", 20);
    final ArrayList<AlignmentMatch> ml = new ArrayList<>();
    ml.add(mA);
    ml.add(mTT);
    ml.add(mGG);
    ml.add(mCC);

    final VariantParams vp = HypothesesComplexTest.getVariantParams(0.5, 0.5, 0.1);

    final ComplexTemplate cot = new ComplexTemplate(DNA.stringDNAtoByte("GGATCGGGGG"), "", 5, 5);
    final HypothesesComplex hyp = HypothesesComplex.makeComplexHypotheses(cot, ml, LogPossibility.SINGLETON, true, vp, null);
    final int readHyp = new EvidenceComplex(hyp, mCC, cot, vp, getChooser()).read();
    assertEquals(2, readHyp);
    final int readHyp2 = new EvidenceComplex(hyp, mTT, cot, vp, getChooser()).read();
    assertEquals(4, readHyp2);
  }

  public void testErrorZeroWhenAllMatchIsSameAsReference() throws Exception {
    final AlignmentMatch mG = HypothesesComplexTest.match("G", 20);
    final ArrayList<AlignmentMatch> ml = new ArrayList<>();
    for (int i = 0; i < 10; i++) {
      ml.add(mG);
    }
    final VariantParams vp = HypothesesComplexTest.getVariantParams(0.5, 0.5, 0.1);
    final ComplexTemplate cot = new ComplexTemplate(DNA.stringDNAtoByte("GGATCGGGGG"), "", 5, 6);
    final HypothesesComplex hyp = HypothesesComplex.makeComplexHypotheses(cot, ml, LogPossibility.SINGLETON, true, vp, null);
    final EvidenceComplex dc = new EvidenceComplex(hyp, mG, cot, vp, getChooser());
    assertEquals(0.0, dc.error());
    assertTrue(dc.integrity());
  }

  public void testErrorRegression() throws Exception {
    final AlignmentMatch mG = HypothesesComplexTest.match("TTTTT", 20);
    final ArrayList<AlignmentMatch> ml = new ArrayList<>();
    ml.add(mG);
    final VariantParams vp = HypothesesComplexTest.getVariantParams(0.5, 0.5, 0.1);
    final ComplexTemplate cot = new ComplexTemplate(DNA.stringDNAtoByte("GGATCGGGGG"), "", 1, 6);
    final HypothesesComplex hyp = HypothesesComplex.makeComplexHypotheses(cot, ml, LogPossibility.SINGLETON, true, vp, null);
    final EvidenceComplex dc = new EvidenceComplex(hyp, mG, cot, vp, getChooser());
    assertEquals(6.498e-6, Math.exp(dc.sumLn()), 1e-6);
    assertEquals(0.99999476, dc.error(), 1e-6);
    assertTrue(dc.integrity());
  }

  public void testReadHypNull() throws Exception {
    final AlignmentMatch mG = HypothesesComplexTest.match("TTTTT", 20);
    final AlignmentMatch mG2 = HypothesesComplexTest.match("CCCC", 20);
    final ArrayList<AlignmentMatch> ml = new ArrayList<>();
    ml.add(mG);
    final VariantParams vp = HypothesesComplexTest.getVariantParams(0.5, 0.5, 0.1);
    final ComplexTemplate cot = new ComplexTemplate(DNA.stringDNAtoByte("GGATCGGGGG"), "", 1, 6);
    final HypothesesComplex hyp = HypothesesComplex.makeComplexHypotheses(cot, ml, LogPossibility.SINGLETON, true, vp, null);
    final EvidenceComplex dc = new EvidenceComplex(hyp, mG2, cot, vp, getChooser());
    assertEquals(EvidenceInterface.NOT_A_HYPOTHESIS, dc.read());
  }

  public void testObviousSingleNtInsertionMatchesHeteroCall() throws Exception {
    final AlignmentMatch mA = HypothesesComplexTest.match("A", 20);
    final AlignmentMatch mTT = HypothesesComplexTest.match("TT", 20);
    final ArrayList<AlignmentMatch> ml = new ArrayList<>();
    for (int k = 0; k < 20; k++) {
      ml.add(mA);
      ml.add(mTT);
    }
    final VariantParams vp = HypothesesComplexTest.getVariantParams(0.5, 0.5, 0.1);

    final ComplexTemplate cot = new ComplexTemplate(DNA.stringDNAtoByte("GGATCGGGGG"), "", 5, 5);
    final HypothesesComplex hyp = HypothesesComplex.makeComplexHypotheses(cot, ml, LogPossibility.SINGLETON, true, vp, null);
    //System.err.println(hyp);
    final StatisticsComplex statistics = new StatisticsComplex(hyp.description(), cot.getLength());
    final ModelInterface<?> m = getModel(hyp, statistics);
    double sum = 0.0;
    for (final AlignmentMatch match : ml) {
      final EvidenceComplex evidence = new EvidenceComplex(hyp, match, cot, vp, getChooser());
      sum += evidence.sumLn();
      m.increment(evidence);
    }
    assertEquals(-469.146, sum, 0.001);
    checkPosterior(m, new Hyp("", -124.665),
        new Hyp("A", -110.021),
        new Hyp("TT", -93.5166),
        new Hyp(":A", -709.705 - sum),
        new Hyp("A:TT", -643.345 - sum),
        new Hyp(":TT", -682.155 - sum));
  }

  public void testBugLongHyp() throws Exception {
    final AlignmentMatch mA = HypothesesComplexTest.match("A", 20);
    final AlignmentMatch mTT = HypothesesComplexTest.match("TTTTTTTTTT", 20);
    final ArrayList<AlignmentMatch> ml = new ArrayList<>();
    ml.add(mTT);
    for (int k = 0; k < 20; k++) {
      ml.add(mA);
    }
    final VariantParams vp = HypothesesComplexTest.getVariantParams(0.5, 0.5, 0.1);

    final ComplexTemplate cot = new ComplexTemplate(DNA.stringDNAtoByte("GGATCGGGGG"), "", 5, 5);
    final HypothesesComplex hyp = HypothesesComplex.makeComplexHypotheses(cot, ml, LogPossibility.SINGLETON, true, vp, null);
    //System.err.println(hyp);
    final StatisticsComplex statistics = new StatisticsComplex(hyp.description(), cot.getLength());
    final ModelInterface<?> m = getModel(hyp, statistics);
    double sum = 0.0;
    for (final AlignmentMatch match : ml) {
      final EvidenceComplex evidence = new EvidenceComplex(hyp, match, cot, vp, getChooser());
      sum += evidence.sumLn();
      m.increment(evidence);
    }
    assertEquals(-223.012, sum, 0.001);
    checkPosterior(m, new Hyp("", -124.665),
        new Hyp("A", -110.021),
        new Hyp("TTTTTTTTTT", -93.5166),
        new Hyp(":A", -709.705 - sum),
        new Hyp("A:TT", -643.345 - sum),
        new Hyp(":TT", -682.155 - sum));
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
    return HypothesesComplexTest.match(sb.toString(), 20);
  }

  public void testRandom() throws Exception {
    for (int j = 0; j < 10; j++) {
      final VariantParams vp = HypothesesComplexTest.getVariantParams(RANDOM.nextDouble() / 10.0, RANDOM.nextDouble() / 10.0, RANDOM.nextDouble() / 10.0);
      final ArrayList<AlignmentMatch> ml = new ArrayList<>();
      for (int k = 0; k < RANDOM.nextInt(500); k++) {
        final AlignmentMatch rm = randomMatch();
        ml.add(rm);
      }

      final HypothesesComplex hyp = HypothesesComplex.makeComplexHypotheses(HypothesesComplexTest.COMPLEX_TEMPLATE, ml, LogPossibility.SINGLETON, true, vp, null);
      final StatisticsComplex statistics = new StatisticsComplex(hyp.description(), HypothesesComplexTest.COMPLEX_TEMPLATE.getLength());
      final ModelInterface<?> m = getModel(hyp, statistics);
      for (final AlignmentMatch match : ml) {
        m.increment(new EvidenceComplex(hyp, match, HypothesesComplexTest.COMPLEX_TEMPLATE, vp, getChooser()));
      }
      m.freeze();
      assertEquals(0, statistics.placedUnmappedCount());
      for (int i = 0 ; i < m.size(); i++) {
        assertFalse(Double.isInfinite(m.posteriorLn0(i)));
        assertFalse(Double.isNaN(m.posteriorLn0(i)));
      }
    }
  }

  // Test probability = 0.0 case
  public void testFromConvergence() throws InvalidParamsException, IOException {
    final ArrayList<AlignmentMatch> ml = new ArrayList<>();

    for (int i = 0; i < 8; i++) {
      ml.add(HypothesesComplexTest.matchi(20));
    }
    final AlignmentMatch cc = HypothesesComplexTest.match("AA", 20);
    for (int i = 0; i < 15; i++) {
      ml.add(cc);
    }
    ml.add(HypothesesComplexTest.match("TGCTAA", 20));
    //System.err.println(ml);
    final GenomePriorParams vpb = new GenomePriorParamsBuilder().create();
    final VariantParams vp = new VariantParamsBuilder().genomePriors(vpb).defaultQuality(20).maxAmbiguity(100.0).maxCoverageFilter(new StaticThreshold(100)).create();
    //System.err.println(ip.toString());
    final HypothesesComplex hyp = HypothesesComplex.makeComplexHypotheses(HypothesesComplexTest.COMPLEX_TEMPLATE, ml, LogPossibility.SINGLETON, true, vp, null);
    final StatisticsComplex statistics = new StatisticsComplex(hyp.description(), HypothesesComplexTest.COMPLEX_TEMPLATE.getLength());
    double sum = 0.0;
    final ModelInterface<?> m = getModel(hyp, statistics);
    for (final AlignmentMatch match : ml) {
      final EvidenceComplex evidence = new EvidenceComplex(hyp, match, HypothesesComplexTest.COMPLEX_TEMPLATE, vp, getChooser());
      sum += evidence.sumLn();
      m.increment(evidence);
    }
    assertEquals(-9.726, sum, 0.001);
    checkPosterior(m, new Hyp("", -144.332),
        new Hyp("AA", -45.1496),
        new Hyp("TGCTAA", -165.2895),
        new Hyp(":AA", -120.100 - sum),
        new Hyp("AA:TGCTAA", -168.462 - sum),
        new Hyp(":TGCTAA", -230.298 - sum));
  }

  // Test probability = 0.0 case including non-identity posterior
  public void testFromConvergence2() throws InvalidParamsException, IOException {
    final ArrayList<AlignmentMatch> ml = new ArrayList<>();

    for (int i = 0; i < 8; i++) {
      ml.add(HypothesesComplexTest.matchi(20));
    }
    final AlignmentMatch cc = HypothesesComplexTest.match("AA", 20);
    for (int i = 0; i < 15; i++) {
      ml.add(cc);
    }
    ml.add(HypothesesComplexTest.match("TGCTAA", 20));

    final GenomePriorParams vpb = new GenomePriorParamsBuilder().create();
    final VariantParams vp = new VariantParamsBuilder().genomePriors(vpb).defaultQuality(20).create();
    //System.err.println(ip.toString());
    final HypothesesComplex hyp = HypothesesComplex.makeComplexHypotheses(HypothesesComplexTest.COMPLEX_TEMPLATE, ml, LogPossibility.SINGLETON, true, vp, null);
    final StatisticsComplex statistics = new StatisticsComplex(hyp.description(), HypothesesComplexTest.COMPLEX_TEMPLATE.getLength());
    final ModelInterface<?> m = getModel(hyp, statistics);
    double sum = 0.0;
    for (final AlignmentMatch match : ml) {
      final EvidenceComplex evidence = new EvidenceComplex(hyp, match, HypothesesComplexTest.COMPLEX_TEMPLATE, vp, getChooser());
      sum += evidence.sumLn();
      m.increment(evidence);
    }
    assertEquals(-9.726, sum, 0.001);
    checkPosterior(m, new Hyp("", -144.332),
        new Hyp("AA", -45.1496),
        new Hyp("TGCTAA", -165.2895),
        new Hyp(":AA", -120.100 - sum),
        new Hyp("AA:TGCTAA", -168.462 - sum),
        new Hyp(":TGCTAA", -230.298 - sum));
  }

  // Test case that gave Infinite posterior
  public void testBug() throws InvalidParamsException, IOException {
    final ArrayList<AlignmentMatch> ml = new ArrayList<>();

    for (int i = 0; i < 9; i++) {
      ml.add(HypothesesComplexTest.matchi(20));
    }
    ml.add(HypothesesComplexTest.match("TGATA", 20));

    final GenomePriorParams vpb = new GenomePriorParamsBuilder().create();
    final VariantParams vp = new VariantParamsBuilder().callLevel(VariantOutputLevel.ALL).genomePriors(vpb).defaultQuality(20).create();
    //System.err.println(ip.toString());
    final HypothesesComplex hyp = HypothesesComplex.makeComplexHypotheses(HypothesesComplexTest.COMPLEX_TEMPLATE, ml, LogPossibility.SINGLETON, true, vp, null);
    final StatisticsComplex statistics = new StatisticsComplex(hyp.description(), HypothesesComplexTest.COMPLEX_TEMPLATE.getLength());
    final ModelInterface<?> m = getModel(hyp, statistics);
    double sum = 0.0;
    for (final AlignmentMatch match : ml) {
      final EvidenceComplex evidence = new EvidenceComplex(hyp, match, HypothesesComplexTest.COMPLEX_TEMPLATE, vp, getChooser());
      sum += evidence.sumLn();
      m.increment(evidence);
    }
    assertEquals(-3.660, sum, 0.001);
    checkPosterior(m, new Hyp("", -15.1521),
        new Hyp("TGATA", -40.358),
        new Hyp(":TGATA", -59.325 - sum)
        );
  }

  // systematically check priors for varying lengths for insertions
  public void testPriorsBig() {
    final ArrayList<AlignmentMatch> ml = new ArrayList<>();
    ml.add(HypothesesComplexTest.match("", 20));
    ml.add(HypothesesComplexTest.match("A", 20));
    ml.add(HypothesesComplexTest.match("AA", 20));
    ml.add(HypothesesComplexTest.match("AAA", 20));
    ml.add(HypothesesComplexTest.match("AAAA", 20));
    ml.add(HypothesesComplexTest.match("AAAAA", 20));
    ml.add(HypothesesComplexTest.match("AAAAAA", 20));
    ml.add(HypothesesComplexTest.match("AAAAAAA", 20));

    final VariantParams vp = VariantParams.builder().genomePriors(GenomePriorParams.builder().create()).create();
    //System.err.println(ip.toString());
    final HypothesesComplex hyp = HypothesesComplex.makeComplexHypotheses(HypothesesComplexTest.COMPLEX_TEMPLATE, ml, LogPossibility.SINGLETON, true, vp, null);
    final StatisticsComplex statistics = new StatisticsComplex(hyp.description(), HypothesesComplexTest.COMPLEX_TEMPLATE.getLength());
    final ModelInterface<?> m = getModel(hyp, statistics);
    double sum = 0.0;
    for (final AlignmentMatch match : ml) {
      final EvidenceComplex evidence = new EvidenceComplex(hyp, match, HypothesesComplexTest.COMPLEX_TEMPLATE, vp, getChooser());
      sum += evidence.sumLn();
      m.increment(evidence);
    }
    assertEquals(-3.523, sum, 0.001);

    checkPosterior(m, new Hyp("", -90.530),
        new Hyp("A", -82.300),
        new Hyp("AA", -70.678),
        new Hyp("AAA", -63.549),
        new Hyp("AAAA", -58.358),
        new Hyp("AAAAA", -58.296),
        new Hyp("AAAAAA", -60.071),
        new Hyp("AAAAAAA", -66.942),
        new Hyp(":A", -110.271 - sum),
        new Hyp("A:AA", -100.103 - sum),
        new Hyp("AA:AAA", -95.490 - sum),
        new Hyp("AAA:AAAA", -92.935 - sum),
        new Hyp("AAAA:AAAAA", -93.842 - sum),
        new Hyp("AAAAA:AAAAAA", -98.377 - sum),
        new Hyp("AAAAAA:AAAAAAA", -106.742 - sum),
        new Hyp(":AA", -97.827 - sum),
        new Hyp("A:AAA", -91.086 - sum),
        new Hyp("AA:AAAA", -87.780 - sum),
        new Hyp("AAA:AAAAA", -89.691 - sum),
        new Hyp("AAAA:AAAAAA", -91.773 - sum),
        new Hyp("AAAAA:AAAAAAA", -100.807 - sum),
        new Hyp(":AAA", -89.835 - sum),
        new Hyp("A:AAAA", -84.161 - sum),
        new Hyp("AA:AAAAA", -85.173 - sum),
        new Hyp("AAA:AAAAAA", -88.150 - sum),
        new Hyp("AAAA:AAAAAAA", -94.640 - sum),
        new Hyp(":AAAA", -83.799 - sum),
        new Hyp("A:AAAAA", -82.334 - sum),
        new Hyp("AA:AAAAAA", -84.333 - sum),
        new Hyp("AAA:AAAAAAA", -91.651 - sum),
        new Hyp(":AAAAA", -82.802 - sum),
        new Hyp("A:AAAAAA", -82.270 - sum),
        new Hyp("AA:AAAAAAA", -88.569 - sum),
        new Hyp(":AAAAAA", -83.594 - sum),
        new Hyp("A:AAAAAAA", -87.311 - sum),
        new Hyp(":AAAAAAA", -89.464));
  }

  public void testSomeCGStuff() throws Exception {
    final SAMRecord sam = new SAMRecord(new SAMFileHeader());
    sam.setAlignmentStart(1);
    sam.setReadString("tatgaaattctgggttgaaaatttttaagaatg");
    sam.setCigarString("21=2I4N10=");
    sam.setAttribute(SamUtils.CG_SUPER_CIGAR, "5=2B18=2I4N10=");
    sam.setAttribute(SamUtils.CG_READ_DELTA, "AA");
    sam.setAttribute(SamUtils.CG_OVERLAP_QUALITY, "55");
    sam.setFlags(73);

    final AlignmentMatch ma = new AlignmentMatch(new VariantAlignmentRecord(sam), "AA", "55", 0, 0, 2, 20);

    final ArrayList<AlignmentMatch> ml = new ArrayList<>();
    ml.add(ma);
    final VariantParams vp = new VariantParamsBuilder().genomePriors(new GenomePriorParamsBuilder().create()).maxCoverageFilter(new StaticThreshold(10)).create();
    final ComplexTemplate cot = new ComplexTemplate(DNA.stringDNAtoByte("gatatgaaattctgggttgaaaattcttttaagaatgttgaatattggcccccacggtcttctggcttgtagggtttctgcagag"), "", 5, 7);

    final DefaultMachineErrorChooser mec = new DefaultMachineErrorChooser("cg_test_errors-080412");
    final HypothesesComplex hyp = HypothesesComplex.makeComplexHypotheses(cot, ml, LogPossibility.SINGLETON, false, vp, null);
    final StatisticsComplex statistics = new StatisticsComplex(hyp.description(), cot.getLength());
    final ModelInterface<?> m = getModel(hyp, statistics);
    double sum = 0.0;
    for (final AlignmentMatch match : ml) {
      final EvidenceComplex evidence = new EvidenceComplex(hyp, match, cot, vp, mec);
      sum += evidence.sumLn();
      m.increment(evidence);
    }
    assertEquals(-15.29198, sum, 0.0001);
    checkPosterior(m,  new Hyp("AA:AA", -9.773),
        new Hyp("GA:GA", -0.015989),
        new Hyp("AA:GA", -17.304));

    final SAMRecord sam3 = new SAMRecord(new SAMFileHeader());
    sam3.setAlignmentStart(1);
    sam3.setReadString("tatgaaattctgggttgaaaatttttaagaatg");
    sam3.setCigarString("21=2I4N10=");
    sam3.setAttribute(SamUtils.CG_SUPER_CIGAR, "5=2B18=2I4N10=");
    sam3.setAttribute(SamUtils.CG_READ_DELTA, "AA");
    sam3.setAttribute(SamUtils.CG_OVERLAP_QUALITY, "55");
    sam3.setFlags(137);

    final AlignmentMatch ma3 = new AlignmentMatch(new VariantAlignmentRecord(sam3), "AA", "55", 0, 0, 2, 20);

    final ArrayList<AlignmentMatch> ml3 = new ArrayList<>();
    ml3.add(ma3);
    final VariantParams vp3 = new VariantParamsBuilder().genomePriors(new GenomePriorParamsBuilder().create()).maxCoverageFilter(new StaticThreshold(10)).create();
    final ComplexTemplate cot3 = new ComplexTemplate(DNA.stringDNAtoByte("gatatgaaattctgggttgaaaattcttttaagaatgttgaatattggcccccacggtcttctggcttgtagggtttctgcagag"), "", 5, 7);

    final DefaultMachineErrorChooser mec3 = new DefaultMachineErrorChooser("cg_test_errors-080412");
    final HypothesesComplex hyp3 = HypothesesComplex.makeComplexHypotheses(cot3, ml3, LogPossibility.SINGLETON, false, vp3, null);
    final StatisticsComplex statistics3 = new StatisticsComplex(hyp3.description(), cot3.getLength());
    final ModelInterface<?> m3 = getModel(hyp3, statistics3);
    double sum1 = 0.0;
    for (final AlignmentMatch match : ml3) {
      final EvidenceComplex evidence = new EvidenceComplex(hyp3, match, cot3, vp3, mec3);
      sum1 += evidence.sumLn();
      m3.increment(evidence);
    }
    assertEquals(-53.9506, sum1, 0.0001);
    checkPosterior(m3, new Hyp("AA:AA", -8.744),
        new Hyp("GA:GA", -0.06448),
        new Hyp("AA:GA", -17.305));

  }

  public void testSoftClippingAtStart() {
    final SAMFileHeader sfh = new SAMFileHeader();
    //                     SSSSSSSSSSSSSSSSSSSSSSSSSSS
    final String readNt = "TTTTTTTTTTGGAGGTGTCCCCAAAATGCCATCAGATGTCATTCACACAATGTATATCTGCACATTATTCCAATACAAGGCAAAGGGGTCTCACATCTGTTAACCGAGTA";
    final SAMRecord sam = new SAMRecord(sfh);
    sam.setAlignmentStart(137);
    sam.setReadString(readNt);
    sam.setCigarString("27S83=");
    sam.setFlags(0);

    final AlignmentMatch ma = new AlignmentMatch(new VariantAlignmentRecord(sam), readNt, null, 0, 0, readNt.length(), 0);

    final String t = "TCATGAATCAGAATCTCATCTTGTAATATTCCAGTGTTCTGTTCTTCAGAAAGTTGCTTCTGAGTATCATCTTGTTCATCACTAGAAAAAAAATTAATTTTCATGAAATACTGGAGGTGTCCCCAAAATG"
                + "ATCTGCG" //deleted region. Ambiguous as to whether starting G is deleted, or the G at the end of this region.
                + "CCATCAGATGTCATTCACACAATGTATATCTGCACATTATTCCAATACAAGGCAAAGGGGTCTCACATCTGTTAACCGAGTATCCCCAACCATGCTGGCACCAGGGACTGGTTTTGCGGAAGATAATTTTTCCAGGAGCCTGAG"; //the start of this region is where the read SHOULD align.

    final ComplexTemplate cot = new ComplexTemplate(DNA.stringDNAtoByte(t), "chrYO", 130, 137);

    final VariantParams vp = new VariantParamsBuilder().genomePriors(new GenomePriorParamsBuilder().create()).maxCoverageFilter(new StaticThreshold(10)).create();
    final MachineErrorChooserInterface mec = getChooser();

    final Match delmatch = new AlleleAsReadMatch(new byte[0]);
    final Match refmatch = new AlleleAsReadMatch(DnaUtils.encodeString("ATCTGCG"));
    final List<Match> durmatches = new ArrayList<>();
    durmatches.add(delmatch);
    durmatches.add(refmatch);

    final DescriptionComplex d = new DescriptionComplex(durmatches) {
    };
    final PossibilityArithmetic arith = LogPossibility.SINGLETON;
    final MockHypotheses<DescriptionComplex> mh = new MockHypotheses<>(d, LogPossibility.SINGLETON, false, new double[] {arith.poss2Prob(-12.132), arith.poss2Prob(-0.279), arith.poss2Prob(-13.725)}, 1);

    final EvidenceComplex evidence = new EvidenceComplex(mh, ma, cot, vp, mec);
    final Double delPoss = arith.prob2Poss(evidence.probability(0));
    final Double refPoss = arith.prob2Poss(evidence.probability(1));

    assertTrue("" + delPoss, delPoss > -1);   //when the start positions are wrong, this value is ~-14, which is wrong. Should be almost 0.
    assertTrue(delPoss > refPoss); //more negative is less probable, we want the del prob to be most likely.
  }

  private AlignmentMatch unmappedMatch() {
    final SAMRecord sam = new SAMRecord(new SAMFileHeader());
    sam.setAlignmentStart(1);
    sam.setReadString("CCCCCGGGGG");
    sam.setReadUnmappedFlag(true);
    return new AlignmentMatch(new VariantAlignmentRecord(sam), "", null, RANDOM.nextInt(64), 0, 0, 10);
  }

  public void testUnmapped() throws Exception {
    final VariantParams vp = HypothesesComplexTest.getVariantParams(RANDOM.nextDouble() / 10.0, RANDOM.nextDouble() / 10.0, RANDOM.nextDouble() / 10.0);
    final AlignmentMatch match = unmappedMatch();
    final HypothesesComplex hyp = HypothesesComplex.makeComplexHypotheses(HypothesesComplexTest.COMPLEX_TEMPLATE, Collections.singletonList(match), LogPossibility.SINGLETON, true, vp, null);
    final StatisticsComplex statistics = new StatisticsComplex(hyp.description(), HypothesesComplexTest.COMPLEX_TEMPLATE.getLength());
    final ModelInterface<?> m = getModel(hyp, statistics);
    m.increment(new EvidenceComplex(hyp, match, HypothesesComplexTest.COMPLEX_TEMPLATE, vp, getChooser()));
    m.freeze();
    for (int i = 0 ; i < m.size(); i++) {
      assertFalse(Double.isInfinite(m.posteriorLn0(i)));
      assertFalse(Double.isNaN(m.posteriorLn0(i)));
    }
    assertEquals(1, statistics.placedUnmappedCount());
  }

}
