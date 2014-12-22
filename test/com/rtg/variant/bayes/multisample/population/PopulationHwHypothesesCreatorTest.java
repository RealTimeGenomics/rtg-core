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



import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import com.rtg.mode.DNA;
import com.rtg.sam.SamRangeUtils;
import com.rtg.tabix.TabixIndexer;
import com.rtg.tabix.UnindexableDataException;
import com.rtg.util.InvalidParamsException;
import com.rtg.util.intervals.RegionRestriction;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.FileUtils;
import com.rtg.util.test.BgzipFileHelper;
import com.rtg.util.test.FileHelper;
import com.rtg.util.test.NanoRegression;
import com.rtg.variant.AlleleCountsFileConverter;
import com.rtg.variant.GenomePriorParams;
import com.rtg.variant.GenomePriorParamsBuilder;
import com.rtg.variant.bayes.complex.ComplexTemplate;
import com.rtg.variant.bayes.complex.DescriptionComplex;
import com.rtg.variant.bayes.complex.HypothesesComplex;
import com.rtg.variant.bayes.complex.HypothesesComplexTest;
import com.rtg.variant.bayes.multisample.HaploidDiploidHypotheses;
import com.rtg.variant.bayes.snp.HypothesesNone;
import com.rtg.variant.bayes.snp.HypothesesPrior;
import com.rtg.variant.bayes.snp.ModelSnpFactory;
import com.rtg.variant.match.AlignmentMatch;

import junit.framework.TestCase;

/**
 */
public class PopulationHwHypothesesCreatorTest extends TestCase {

  private static final String TAB = "\t";
  private static final String VCF = ""
      + "##fileformat=VCFv4.0" + "\n"
      + "##source=BCM:SNPTools:hapfuse" + "\n"
      + "##reference=1000Genomes-NCBI37" + "\n"
      + "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" + "\n"
      + "##FORMAT=<ID=AP,Number=2,Type=Float,Description=\"Allelic Probability, P(Allele=1|Haplotype)\">" + "\n"
      + "#CHROM" + TAB + "POS" + TAB + "ID" + TAB + "REF" + TAB + "ALT" + TAB + "QUAL" + TAB + "FILTER" + TAB + "INFO" + TAB + "FORMAT" + TAB + "HG00096" + TAB + "HG00097" + TAB + "HG00099" + TAB + "HG00100" + TAB + "HG00101" + "\n"
      + "1" + TAB + "10" + TAB + "rs117577454" + TAB + "C" + TAB + "G" + TAB + "100" + TAB + "PASS" + TAB + "." + TAB + "GT:AP" + TAB + "0|0:0.075,0.060" + TAB + "1|0:0.640,0.000" + TAB + "0|0:0.115,0.015" + TAB + "0|0:0.020,0.000" + TAB + "0|0:0.060,0.000" + "\n"
      + "1" + TAB + "12" + TAB + "rs117577454" + TAB + "CT" + TAB + "GT" + TAB + "100" + TAB + "PASS" + TAB + "." + TAB + "GT:AP" + TAB + "0|0:0.075,0.060" + TAB + "1|0:0.640,0.000" + TAB + "0|0:0.115,0.015" + TAB + "0|0:0.020,0.000" + TAB + "0|0:0.060,0.000" + "\n"
      + "1" + TAB + "15" + TAB + "rs117577454" + TAB + "." + TAB + "GT" + TAB + "100" + TAB + "PASS" + TAB + "." + TAB + "GT:AP" + TAB + "0|0:0.075,0.060" + TAB + "1|0:0.640,0.000" + TAB + "0|0:0.115,0.015" + TAB + "0|0:0.020,0.000" + TAB + "0|0:0.060,0.000" + "\n"
      + "1" + TAB + "18" + TAB + "rs117577454" + TAB + "C" + TAB + "G" + TAB + "100" + TAB + "BLAH" + TAB + "." + TAB + "GT:AP" + TAB + "0|1:0.075,0.060" + TAB + "1|0:0.640,0.000" + TAB + "0|0:0.115,0.015" + TAB + "0|0:0.020,0.000" + TAB + "0|0:0.060,0.000" + "\n"
      + "1" + TAB + "19" + TAB + "rs117577454" + TAB + "CTGAGG" + TAB + "C,CTGAAG" + TAB + "100" + TAB + "PASS" + TAB + "XRX" + TAB + "GT:AP" + TAB + "0|0:0.075,0.060" + TAB + "1|2:0.640,0.000" + TAB + "0|0:0.115,0.015" + TAB + "0|0:0.020,0.000" + TAB + "0|0:0.060,0.000" + "\n";

  NanoRegression mNano = null;

  @Override
  public void setUp() {
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

  public void test() throws IOException, UnindexableDataException, InvalidParamsException {
    final File dir = FileUtils.createTempDir("populationhw", "test");
    try {
      final File out = new File(dir, "snps.vcf.gz");
      FileHelper.stringToGzFile(VCF, out);

      final File alleleCountFile = new File(dir, "allele.ac");

      final AlleleCountsFileConverter blah = new AlleleCountsFileConverter();
      blah.convert(out, alleleCountFile);

      new TabixIndexer(out, new File(dir, "snps.vcf.gz.tbi")).saveVcfIndex();
      final GenomePriorParams genomePriors = GenomePriorParams.builder().create();
      final ModelSnpFactory haploid = new ModelSnpFactory(genomePriors, true);
      final ModelSnpFactory diploid = new ModelSnpFactory(genomePriors, false);
      final PopulationHwHypothesesCreator pwh = new PopulationHwHypothesesCreator(alleleCountFile, haploid, diploid, null);
      final HaploidDiploidHypotheses<?> hypotheses = pwh.getSnpHypotheses("1", 9);

      assertNotNull(hypotheses);
      assertEquals(true, hypotheses.haploid().haploid());
      assertEquals(false, hypotheses.diploid().haploid());

      assertNull(pwh.getSnpHypotheses("1", 11));
      assertNull(pwh.getSnpHypotheses("1", 14));

      mNano.check("pophwhyptest", hypotheses.haploid().toString() + hypotheses.diploid().toString());
      /*
      System.err.println("1 " + pwh.getCounts("1", 10, 12));
      System.err.println("2 " + pwh.getCounts("1", 1, 9));
      System.err.println("3 " + pwh.getCounts("1", 20, 30));
      System.err.println("4 " + pwh.getCounts("2", 10, 12));
      System.err.println("5 " + pwh.getCounts("1", 15, 30));
      System.err.println("6 " + pwh.getCounts("1", 1, 10));
       */
    } finally {
      FileHelper.deleteAll(dir);
    }
  }

  public void testGetCounts() throws IOException, UnindexableDataException, InvalidParamsException {
    final File dir = FileUtils.createTempDir("populationhw", "test");
    try {
      final File out = new File(dir, "snps.vcf.gz");
      FileHelper.stringToGzFile(VCF, out);

      final File alleleCountFile = new File(dir, "allele.ac");

      final AlleleCountsFileConverter blah = new AlleleCountsFileConverter();
      blah.convert(out, alleleCountFile);

      new TabixIndexer(out, new File(dir, "snps.vcf.gz.tbi")).saveVcfIndex();
      final GenomePriorParams genomePriors = GenomePriorParams.builder().create();
      final ModelSnpFactory haploid = new ModelSnpFactory(genomePriors, true);
      final ModelSnpFactory diploid = new ModelSnpFactory(genomePriors, false);
      final PopulationHwHypothesesCreator pwh = new PopulationHwHypothesesCreator(alleleCountFile, haploid, diploid, null);
      assertEquals(0, pwh.getCounts("1", 0, 5).size());
      assertEquals(0, pwh.getCounts("1", 0, 9).size());
      assertEquals(1, pwh.getCounts("1", 0, 10).size());
      assertEquals(2, pwh.getCounts("1", 0, 15).size());        //this changed from 3 to 2 due to disallowing vcf lines with a reference of '.'
      assertEquals(2, pwh.getCounts("1", 0, 14).size());
      assertEquals(1, pwh.getCounts("1", 16, 19).size());
      assertEquals(1, pwh.getCounts("1", 18, 30).size());
      assertEquals(1, pwh.getCounts("1", 23, 30).size());
      assertEquals(0, pwh.getCounts("1", 24, 30).size());
      assertEquals(0, pwh.getCounts("1", 25, 30).size());
    } finally {
      FileHelper.deleteAll(dir);
    }
  }

  // TODO Following test should be enabled after fixing TODO in HypothesesComplex.createDescription method
  //  public void testGetComplexHypothesesPruning() {
  //    final String ref = "CCCACTGTTGC";
  //
  //    final ArrayList<AlignmentMatch> ml = new ArrayList<>();
  //    for (int i = 0; i <= 934; i++) {
  //      ml.add(HypothesesComplexTest.match(ref, 40));
  //    }
  //    for (int i = 0; i <= 6; i++) {
  //      ml.add(HypothesesComplexTest.match("CCACTGTTGC", 40));
  //    }
  //    for (int i = 0; i <= 5; i++) {
  //      ml.add(HypothesesComplexTest.match("CCCACTGCTGC", 40));
  //    }
  //    for (int i = 0; i <= 5; i++) {
  //      ml.add(HypothesesComplexTest.match("CCCACTGTCGC", 40));
  //    }
  //    for (int i = 0; i <= 5; i++) {
  //      ml.add(HypothesesComplexTest.match("CCCCACTGTTGC", 40));
  //    }
  //    ml.add(HypothesesComplexTest.match("CCCACTAGTTGC", 40));
  //    ml.add(HypothesesComplexTest.match("CCCACTAGTTGC", 40));
  //
  //    ml.add(HypothesesComplexTest.match("CCAACTTGTTGC", 40));
  //    ml.add(HypothesesComplexTest.match("CCACATGTTGC", 40));
  //    ml.add(HypothesesComplexTest.match("CCACTGTGC", 40));
  //    ml.add(HypothesesComplexTest.match("CCATGTTGC", 40));
  //    ml.add(HypothesesComplexTest.match("CCCACATGTTGC", 40));
  //    ml.add(HypothesesComplexTest.match("CCCACCGTTGC", 40));
  //    ml.add(HypothesesComplexTest.match("CCCACTATTGC", 40));
  //    ml.add(HypothesesComplexTest.match("CCCACTGTGC", 40));
  //    ml.add(HypothesesComplexTest.match("CCCACTGTTAC", 40));
  //    ml.add(HypothesesComplexTest.match("CCCACTGTTAGC", 40));
  //    ml.add(HypothesesComplexTest.match("CCCACTGTTGCA", 40));
  //    ml.add(HypothesesComplexTest.match("CCCACTGTTGT", 40));
  //    ml.add(HypothesesComplexTest.match("CCCACTTGTTGC", 40));
  //    ml.add(HypothesesComplexTest.match("CCCACTTGTTGCA", 40));
  //    ml.add(HypothesesComplexTest.match("CCCCCTGTTGC", 40));
  //    ml.add(HypothesesComplexTest.match("CCCGCTGTTGC", 40));
  //    ml.add(HypothesesComplexTest.match("CCCTACTGTTGC", 40));
  //    ml.add(HypothesesComplexTest.match("CCTACTGTTGC", 40));
  //    ml.add(HypothesesComplexTest.match("TCCACTGTTGC", 40));
  //    ml.add(HypothesesComplexTest.match("TCCCACTGTTG", 40));
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
  //        counts.put("CCACATGTTGC", 1100);
  //        counts.put("AAAAAAAAAA", 10);
  //        counts.put("T", 1);
  //        final ArrayList<AlleleCounts> list = new ArrayList<>();
  //        list.add(new AlleleCounts(1, counts, 1, false));
  //        return list;
  //      }
  //    };
  //
  //    final ComplexTemplate cot = new ComplexTemplate(DNA.stringDNAtoByte(ref), "", 0, ref.length());
  //
  //    final HypothesesComplex haploid = HypothesesComplexTest.getSpecifiedScorer(ml, true, ssp, cot, true);
  //    final HypothesesComplex diploid = HypothesesComplexTest.getSpecifiedScorer(ml, false, ssp, cot, true);
  //
  //    assertEquals(haploid.toString(), 10, haploid.description().size());
  //
  //    final HaploidDiploidHypotheses<?> res = PopulationHwHypothesesCreator.getComplexHypotheses(new HaploidDiploidHypotheses<HypothesesComplex>(haploid, diploid));
  //    assertNotNull(res);
  //    assertNotNull(res.haploid());
  //    assertNotNull(res.diploid());
  //    assertTrue(res.haploid() instanceof HypothesesComplex);
  //    assertTrue(res.diploid() instanceof HypothesesComplex);
  //
  //
  //    mNano.check("pophwhypcomplexpruning", haploid.toString() + diploid.toString() + res.haploid().toString() + res.diploid().toString(), true);
  //
  //  }
  //
  //  // This is testing a pathological case that in theory could be improved at some point.
  //  public void testGetComplexHypothesesPruningPathological() {
  //    final String ref = "CCCACTGTTGC";
  //
  //    final ArrayList<AlignmentMatch> ml = new ArrayList<>();
  //    for (int i = 0; i <= 25; i++) {
  //      ml.add(HypothesesComplexTest.match(ref, 40));
  //    }
  //    for (int i = 0; i <= 25; i++) {
  //      ml.add(HypothesesComplexTest.match("CCTACTGATGC", 40));
  //    }
  //
  //    final SiteSpecificPriors ssp = new SiteSpecificPriors() {
  //      @Override
  //      public HaploidDiploidHypotheses<?> getHypotheses(String templateName, int zeroPos) {
  //        return null;
  //      }
  //
  //      @Override
  //      public List<AlleleCounts> getCounts(String templateName, int start, int end) {
  //        HashMap<String, Integer> counts = new HashMap<>();
  //        counts.put("T", 50);
  //        counts.put("C", 50);
  //        final ArrayList<AlleleCounts> list = new ArrayList<>();
  //        list.add(new AlleleCounts(2, counts, 1, false));
  //
  //        counts = new HashMap<>();
  //        counts.put("A", 90);
  //        counts.put("T", 10);
  //
  //        list.add(new AlleleCounts(7, counts, 1, false));
  //        return list;
  //      }
  //    };
  //
  //    final ComplexTemplate cot = new ComplexTemplate(DNA.stringDNAtoByte(ref), "", 0, ref.length());
  //
  //    final HypothesesComplex haploid = HypothesesComplexTest.getSpecifiedScorer(ml, true, ssp, cot, true);
  //    final HypothesesComplex diploid = HypothesesComplexTest.getSpecifiedScorer(ml, false, ssp, cot, true);
  //
  //    assertEquals(haploid.toString(), 4, haploid.description().size());
  //
  //    final HaploidDiploidHypotheses<?> res = PopulationHwHypothesesCreator.getComplexHypotheses(new HaploidDiploidHypotheses<HypothesesComplex>(haploid, diploid));
  //
  //    assertNotNull(res);
  //    assertNotNull(res.haploid());
  //    assertNotNull(res.diploid());
  //    assertTrue(res.haploid() instanceof HypothesesComplex);
  //    assertTrue(res.diploid() instanceof HypothesesComplex);
  //
  //    mNano.check("poppathological", haploid.toString() + diploid.toString() + res.haploid().toString() + res.diploid().toString(), true);
  //
  //  }
  //
  //
  //
  //  public void testGetComplexHypotheses() {
  //    final String ref = "CCCACTGTTGC";
  //
  //    final ArrayList<AlignmentMatch> ml = new ArrayList<>();
  //    for (int i = 0; i <= 934; i++) {
  //      ml.add(HypothesesComplexTest.match(ref, 40));
  //    }
  //    for (int i = 0; i <= 6; i++) {
  //      ml.add(HypothesesComplexTest.match("CCACTGTTGC", 40));
  //    }
  //    for (int i = 0; i <= 5; i++) {
  //      ml.add(HypothesesComplexTest.match("CCCACTGCTGC", 40));
  //    }
  //    for (int i = 0; i <= 5; i++) {
  //      ml.add(HypothesesComplexTest.match("CCCACTGTCGC", 40));
  //    }
  //    for (int i = 0; i <= 5; i++) {
  //      ml.add(HypothesesComplexTest.match("CCCCACTGTTGC", 40));
  //    }
  //    ml.add(HypothesesComplexTest.match("CCCACTAGTTGC", 40));
  //    ml.add(HypothesesComplexTest.match("CCCACTAGTTGC", 40));
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
  //        counts.put("CCACATGTTGC", 1100);
  //        counts.put("AAAAAAAAAA", 10);
  //        counts.put("T", 1);
  //        final ArrayList<AlleleCounts> list = new ArrayList<>();
  //        list.add(new AlleleCounts(1, counts, 1, false));
  //        return list;
  //      }
  //    };
  //
  //    final ComplexTemplate cot = new ComplexTemplate(DNA.stringDNAtoByte(ref), "", 0, ref.length());
  //
  //    final HypothesesComplex haploid = HypothesesComplexTest.getSpecifiedScorer(ml, true, ssp, cot, false);
  //    final HypothesesComplex diploid = HypothesesComplexTest.getSpecifiedScorer(ml, false, ssp, cot, false);
  //
  //    assertEquals(10, haploid.description().size());
  //
  //    final HaploidDiploidHypotheses<?> res = PopulationHwHypothesesCreator.getComplexHypotheses(new HaploidDiploidHypotheses<HypothesesComplex>(haploid, diploid));
  //    assertNotNull(res);
  //    assertNotNull(res.haploid());
  //    assertNotNull(res.diploid());
  //    assertTrue(res.haploid() instanceof HypothesesComplex);
  //    assertTrue(res.diploid() instanceof HypothesesComplex);
  //
  //    mNano.check("pophwhypcomplex", haploid.toString() + diploid.toString() + res.haploid().toString() + res.diploid().toString(), true);
  //  }


  public void testGetComplexHypothesesNovelty() throws Exception {
    final String ref = "CCCACTGTTGC";

    final ArrayList<AlignmentMatch> ml = new ArrayList<>();
    for (int i = 0; i <= 934; i++) {
      ml.add(HypothesesComplexTest.match(ref, 40));
    }
    for (int i = 0; i <= 6; i++) {
      ml.add(HypothesesComplexTest.match("CCACTGTTGC", 40));
    }
    for (int i = 0; i <= 5; i++) {
      ml.add(HypothesesComplexTest.match("CCCACTGCTGC", 40));
    }
    for (int i = 0; i <= 5; i++) {
      ml.add(HypothesesComplexTest.match("CCCACTGTCGC", 40));
    }
    for (int i = 0; i <= 5; i++) {
      ml.add(HypothesesComplexTest.match("CCCCACTGTTGC", 40));
    }
    ml.add(HypothesesComplexTest.match("CCCACTAGTTGC", 40));
    ml.add(HypothesesComplexTest.match("CCCACTAGTTGC", 40));


    final ComplexTemplate cot = new ComplexTemplate(DNA.stringDNAtoByte(ref), "", 0, ref.length());

    final HypothesesComplex haploid = HypothesesComplexTest.getSpecifiedScorerProb(ml, true, null, cot, false);
    final HypothesesComplex diploid = HypothesesComplexTest.getSpecifiedScorerProb(ml, false, null, cot, false);

    assertEquals(6, haploid.description().size());
    final double novelty = 0.02;
    final HaploidDiploidHypotheses<HypothesesPrior<DescriptionComplex>> hyp = new HaploidDiploidHypotheses<>(HypothesesNone.SINGLETON_COMPLEX, haploid, diploid);
    final HaploidDiploidHypotheses<?> res = PopulationHwHypothesesCreator.computeNoveltyPriors(hyp, novelty);
    assertNotNull(res);
    assertNotNull(res.haploid());
    assertNotNull(res.diploid());
    assertEquals(1.0 - novelty, res.haploid().p(res.haploid().reference()));
    assertEquals(1.0 - novelty, res.diploid().p(res.haploid().reference()));
    assertNotNull(res.haploid());
    assertNotNull(res.diploid());

    mNano.check("pophwnovelty", haploid.toString() + diploid.toString() + res.haploid().toString() + res.diploid().toString());
  }

  private static final String VCFHEADER = "##fileformat=VCFv4.1\n"
      + "##contig=<ID=Chr1,length=158337067>\n"
      + "##INFO=<ID=XRX,Number=0,Type=Flag,Description=\"RTG variant was called using complex caller\">\n"
      + "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Combined read depth for variant over all samples\">\n"
      + "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
      + "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">\n"
      + "##FORMAT=<ID=GQ,Number=1,Type=Float,Description=\"Genotype Quality\">\n"
      + "##FORMAT=<ID=RE,Number=1,Type=Float,Description=\"RTG Total Error\">\n"
      + "##FORMAT=<ID=AR,Number=1,Type=Float,Description=\"Ambiguity Ratio\">\n"
      + "##FORMAT=<ID=AB,Number=1,Type=Float,Description=\"Allele Balance\">\n"
      + "##FORMAT=<ID=RQ,Number=1,Type=Float,Description=\"RTG sample quality\">\n"
      + "##FORMAT=<ID=RS,Number=.,Type=String,Description=\"RTG Support Statistics\">\n"
      + "#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  501038  17144783        18033109        99594   99592   99593   664221  96329   363     664195  17074197        89429   17396939        17108826        18143602   372     346     27      18186561        17145716        586     18093560\n".replaceAll(" +", "\t");


  public void testVariantFilterVcf() throws Exception {
    variantFilterCheck(false);
  }
  public void testVariantFilterAlleleCount() throws Exception {
    variantFilterCheck(true);
  }

  private void variantFilterCheck(boolean useAlleleCountFile) throws Exception {
    final String ssps = "Chr1    23996   .       T       C       .       PASS    DP=1782 GT:DP:RE:AR:AB:RQ:GQ:RS 0/0:99:0.900:0.000:1.000:0.0:304.3:T,99,0.900   0/0:32:0.180:0.000:1.000:0.0:103.7:T,32,0.180   0/0:26:0.100:0.000:1.000:0.0:85.8:T,26,0.100    0/0:118:0.507:0.000:1.000:0.0:362.1:T,118,0.507 0/0:85:0.515:0.000:1.000:0.0:262.8:T,85,0.515   0/0:62:0.507:0.000:1.000:0.0:193.5:T,62,0.507   0/0:98:0.447:0.000:1.000:0.0:302.0:T,98,0.447   1/0:103:1.658:0.010:0.883:8.9:8.9:C,12,1.040,T,91,0.618 0/0:95:0.449:0.000:0.989:0.0:263.2:C,1,0.003,T,94,0.446 0/0:98:0.482:0.000:1.000:0.0:302.0:T,98,0.482   0/0:29:0.129:0.000:1.000:0.0:94.8:T,29,0.129    0/0:96:0.398:0.000:1.000:0.0:296.1:T,96,0.398   0/0:21:0.071:0.000:1.000:0.0:70.8:T,21,0.071    0/0:24:0.095:0.000:0.917:0.0:19.7:C,2,0.006,T,22,0.089  1/0:26:0.111:0.000:0.846:25.8:25.8:C,4,0.026,T,22,0.085 0/0:81:0.323:0.000:1.000:0.0:251.0:T,81,0.323   1/0:235:1.339:0.000:0.851:289.0:289.0:C,35,0.161,T,200,1.177    1/0:85:0.473:0.000:0.576:752.2:752.2:C,36,0.262,T,49,0.211      0/0:27:0.352:0.000:0.963:0.0:85.4:G,1,0.013,T,26,0.339  0/0:32:0.147:0.000:1.000:0.0:103.8:T,32,0.147   0/0:90:0.511:0.000:1.000:0.0:277.8:T,90,0.511   0/0:31:0.221:0.000:1.000:0.0:100.6:T,31,0.221\n".replaceAll(" +", "\t")
        + "Chr1    24000   .       A       G       .       PASS    DP=1794;XRX     GT:DP:RE:RQ:GQ  0/0:99:0.743:889.8:889.5        0/1:30:0.148:90.9:90.5  0/0:23:0.116:0.0:71.0   0/1:120:0.645:440.3:440.3       0/1:92:0.510:253.7:253.7        0/0:61:0.463:363.3:8.7  0/0:97:0.459:0.0:294.2  0/0:106:1.664:342.7:342.7       0/0:88:0.392:0.0:261.7  0/1:103:0.478:1222.8:713.0      0/1:32:0.149:103.3:44.4 0/0:95:0.429:361.3:12.5 0/1:20:1.073:10.7:10.7  0/0:25:0.099:46.9:46.9  0/0:26:1.103:0.0:72.9   0/0:87:0.399:313.8:313.8        0/0:230:1.436:0.0:683.8 0/0:92:1.479:1028.0:741.5       0/1:27:0.363:102.3:102.3        0/0:30:0.144:40.2:40.2  0/1:89:3.536:211.2:211.2        0/1:29:0.199:31.7:31.7\n".replaceAll(" +", "\t")
        + "Chr1    24007   .       G       A       .       PASS    DP=1794;XRX     GT:DP:RE:RQ:GQ  1/0:99:0.743:889.8:889.5        0/0:30:0.148:90.9:90.5  0/0:23:0.116:0.0:71.0   0/0:120:0.645:440.3:440.3       0/0:92:0.510:253.7:253.7        1/0:61:0.463:363.3:8.7  0/0:97:0.459:0.0:294.2  1/0:106:1.664:342.7:342.7       0/0:88:0.392:0.0:261.7  0/0:103:0.478:1222.8:713.0      0/0:32:0.149:103.3:44.4 1/0:95:0.429:361.3:12.5 0/0:20:1.073:10.7:10.7  1/0:25:0.099:46.9:46.9  0/0:26:1.103:0.0:72.9   1/0:87:0.399:313.8:313.8        0/0:230:1.436:0.0:683.8 1/0:92:1.479:1028.0:741.5       0/0:27:0.363:102.3:102.3        1/0:30:0.144:40.2:40.2  0/0:89:3.536:211.2:211.2        0/0:29:0.199:31.7:31.7\n".replaceAll(" +", "\t");

    final File tmpDir = FileHelper.createTempDirectory();
    try {
      final File sspf = new File(tmpDir, "ssps.vcf.gz");
      BgzipFileHelper.bytesToBgzipFile((VCFHEADER + ssps).getBytes(), sspf);

      File fileToUse;
      if (useAlleleCountFile) {
        fileToUse = new File(tmpDir, "allele.ac");
        final AlleleCountsFileConverter blah = new AlleleCountsFileConverter();
        blah.convert(sspf, fileToUse);
      } else {
        final TabixIndexer ti = new TabixIndexer(sspf);
        ti.saveVcfIndex();
        fileToUse = sspf;
      }

      final PopulationHwHypothesesCreator pop = new PopulationHwHypothesesCreator(fileToUse, new GenomePriorParamsBuilder().create(), SamRangeUtils.createExplicitReferenceRange(new RegionRestriction("Chr1")));

      final List<AlleleCounts> acs = pop.getCounts("Chr1", 24000, 25000);
      assertEquals(1, acs.size());
      assertEquals(24006, acs.get(0).position());

    } finally {
      FileHelper.deleteAll(tmpDir);
    }
  }

  public void testInsertChunkBoundary() throws Exception {
    final String ssps = ""
        + "Chr1    24000   .       A       GT       .       PASS    DP=1794;XRX     GT:DP:RE:RQ:GQ  0/0:99:0.743:889.8:889.5        0/1:30:0.148:90.9:90.5  0/0:23:0.116:0.0:71.0   0/1:120:0.645:440.3:440.3       0/1:92:0.510:253.7:253.7        0/0:61:0.463:363.3:8.7  0/0:97:0.459:0.0:294.2  0/0:106:1.664:342.7:342.7       0/0:88:0.392:0.0:261.7  0/1:103:0.478:1222.8:713.0      0/1:32:0.149:103.3:44.4 0/0:95:0.429:361.3:12.5 0/1:20:1.073:10.7:10.7  0/0:25:0.099:46.9:46.9  0/0:26:1.103:0.0:72.9   0/0:87:0.399:313.8:313.8        0/0:230:1.436:0.0:683.8 0/0:92:1.479:1028.0:741.5       0/1:27:0.363:102.3:102.3        0/0:30:0.144:40.2:40.2  0/1:89:3.536:211.2:211.2        0/1:29:0.199:31.7:31.7\n".replaceAll(" +", "\t");

    final File tmpDir = FileHelper.createTempDirectory();
    try {
      final File sspf = new File(tmpDir, "ssps.vcf.gz");
      BgzipFileHelper.bytesToBgzipFile((VCFHEADER + ssps).getBytes(), sspf);

      final File alleleCountFile = new File(tmpDir, "allele.ac");

      final AlleleCountsFileConverter blah = new AlleleCountsFileConverter();
      blah.convert(sspf, alleleCountFile);

      final TabixIndexer ti = new TabixIndexer(sspf);
      ti.saveVcfIndex();

      final PopulationHwHypothesesCreator pop = new PopulationHwHypothesesCreator(alleleCountFile, new GenomePriorParamsBuilder().create(), SamRangeUtils.createExplicitReferenceRange(new RegionRestriction("Chr1")));

      List<AlleleCounts> acs = pop.getCounts("Chr1", 24000, 25000);
      assertEquals(0, acs.size());

      acs = pop.getCounts("Chr1", 23000, 24000);
      assertEquals(1, acs.size());      //we want the insert in the bucket before 24000
    } finally {
      FileHelper.deleteAll(tmpDir);
    }
  }
}
