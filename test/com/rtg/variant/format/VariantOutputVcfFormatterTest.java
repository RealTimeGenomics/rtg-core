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
package com.rtg.variant.format;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.StringReader;
import java.util.EnumSet;

import com.rtg.reference.Ploidy;
import com.rtg.relation.GenomeRelationships;
import com.rtg.relation.RelationshipsFileParser;
import com.rtg.util.MathUtils;
import com.rtg.util.cli.CommandLine;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.util.test.NanoRegression;
import com.rtg.variant.StaticThreshold;
import com.rtg.variant.Variant;
import com.rtg.variant.Variant.VariantFilter;
import com.rtg.variant.VariantLocus;
import com.rtg.variant.VariantParams;
import com.rtg.variant.VariantParamsBuilder;
import com.rtg.variant.VariantSample;
import com.rtg.variant.bayes.MockGenotypeMeasure;
import com.rtg.variant.bayes.NoNonIdentityMeasure;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import junit.framework.TestCase;

/**
 */
public class VariantOutputVcfFormatterTest extends TestCase {

  private static final String TAB = "\t";

  private NanoRegression mNano = null;

  @Override
  public void tearDown() throws Exception {
    CommandLine.clearCommandArgs();
    Diagnostic.setLogStream();
    try {
      mNano.finish();
    } finally {
      mNano = null;
    }
  }

  @Override
  public void setUp() {
    Diagnostic.setLogStream();
    mNano = new NanoRegression(this.getClass(), false);
  }

  public static VariantSample createSample(Ploidy diploid, String name, boolean identity, Double posterior, VariantSample.DeNovoStatus deNovo, Double deNovoPosterior) {
    return createSample(diploid, name, identity, posterior, deNovo, deNovoPosterior, null);
  }
  public static VariantSample createSample(Ploidy diploid, String name, boolean identity, Double posterior, VariantSample.DeNovoStatus deNovo, Double deNovoPosterior, Double nonIdentity) {
    return new VariantSample(diploid, name, identity, new NoNonIdentityMeasure(new MockGenotypeMeasure(0, 1, posterior == null ? 0.0 : posterior, nonIdentity == null ? Double.NaN : nonIdentity)), deNovo, deNovoPosterior);
  }

  public void testHeader1() throws IOException {
    //set command line
    try (MemoryPrintStream dps = new MemoryPrintStream()) {
      Diagnostic.setLogStream(dps.printStream());
      CommandLine.setCommandArgs("yo", "--ho", "ho");

      try (MemoryPrintStream mps = new MemoryPrintStream()) {
        final VariantParams params = new VariantParamsBuilder().maxAmbiguity(0.5).regionsFilterBedFile(new File("foo.bed")).create();
        final VariantOutputVcfFormatter formatter = new VariantOutputVcfFormatter(params, "SAMPLE");

        final SAMFileHeader samHeader = new SAMFileHeader();
        SAMReadGroupRecord readGroup = new SAMReadGroupRecord("s");
        readGroup.setPlatform("Illumina");
        readGroup.setSample("the_one_ring");
        samHeader.addReadGroup(readGroup);
        readGroup = new SAMReadGroupRecord("b");
        readGroup.setPlatform("Illumina");
        readGroup.setSample("blah");
        samHeader.addReadGroup(readGroup);
        readGroup = new SAMReadGroupRecord("c");
        readGroup.setPlatform("Illumina");
        readGroup.setSample("rargh");
        samHeader.addReadGroup(readGroup);

        formatter.writeHeader(mps.outputStream(), params, samHeader);

        assertEquals(1, dps.toString().trim().split("\n").length);
        assertTrue(mps.toString(), mps.toString().contains("SAMPLE"));

        mps.reset();

        formatter.writeHeader(mps.outputStream(), params, null);

        final String header = mps.toString().replaceAll("##source=.*", "##source=").replaceAll("##RUN-ID=.*", "##RUN-ID=").replaceAll("##fileDate=20.*", "##fileDate=20");
        mNano.check("vovf-single-header1.vcf", header);

        //chr10   82350   o       C       T       0.8     4       0.400   T       4       0.400
        final VariantSample vs = createSample(Ploidy.DIPLOID, "T:T", false, 0.8 * MathUtils.LOG_10, VariantSample.DeNovoStatus.NOT_DE_NOVO, 0.0);
        vs.setCoverage(4);
        vs.setCoverageCorrection(0.4);

        final VariantLocus locus = new VariantLocus("chr10", 82349, 82350, "C", 'N');
        final Variant v = new Variant(locus, vs);
        v.setNonIdentityPosterior(2.0);

        final String s = formatter.formatCall(v);

        //chr10      82350   .       C       T       8.63       PASS    NS=1;DP=4;AF=1.000;RE=0.400     GT:GQ:RP:DP:RS  1/1:5:0.8:4:T,4,0.400
        assertEquals("chr10\t82350\t.\tC\tT\t9.2\tPASS\t.\tGT:DP:RE:GQ\t1/1:4:0.400:9\n", s);
      }
    }
  }

  public void test1() {
    final VariantParams params = new VariantParamsBuilder().maxAmbiguity(0.5).create();
    final VariantOutputVcfFormatter formatter = new VariantOutputVcfFormatter(params, "SAMPLE");

    //chr10   112109  e       C       C:T     9.5     38      4.622   C       14      1.430   T       24      3.192
    final VariantSample vs = createSample(Ploidy.DIPLOID, "C:T", false, 9.5 * MathUtils.LOG_10, VariantSample.DeNovoStatus.NOT_DE_NOVO, 0.0);
    vs.setCoverage(38);
    vs.setCoverageCorrection(4.622);

    //chr10      112109  .       C       T       41      PASS    NS=1;DP=38;AF=0.632;RE=4.622    GT:GQ:RP:DP:RS  0/1:41:9.5:38:C,14,1.430,T,24,3.192

    final VariantLocus locus = new VariantLocus("chr10", 112108, 112109, "C", 'N');
    final Variant v = new Variant(locus, vs);
    final String s = formatter.formatCall(v);

    assertEquals("chr10\t112109\t.\tC\tT\t.\tPASS\t.\tGT:DP:RE:GQ\t0/1:38:4.622:95\n", s);
  }

  public void test2() {
    final VariantParams params = new VariantParamsBuilder().maxAmbiguity(0.5).create();
    final VariantOutputVcfFormatter formatter = new VariantOutputVcfFormatter(params, "SAMPLE");

    //chr10   115008  e       G       C:G     16.4    81      15.368  C       53      10.713  G       28      4.654
    final VariantSample vs = createSample(Ploidy.DIPLOID, "C:G", false, 16.4 * MathUtils.LOG_10, VariantSample.DeNovoStatus.UNSPECIFIED, null);
    vs.setCoverage(81);
    vs.setCoverageCorrection(15.368);

    //chr10      115008  .       G       C       71      PASS    NS=1;DP=81;AF=0.654;RE=15.368   GT:GQ:RP:DP:RS  1/0:71:16.4:81:C,53,10.713,G,28,4.654
    final VariantLocus locus = new VariantLocus("chr10", 115007, 115008, "G", 'N');
    final Variant v = new Variant(locus, vs);
    final String s = formatter.formatCall(v);

    assertEquals("chr10\t115008\t.\tG\tC\t.\tPASS\t.\tGT:DP:RE:GQ\t1/0:81:15.368:164\n", s);
  }

  public void test3() {
    final VariantParams params = new VariantParamsBuilder().maxAmbiguity(0.5).create();
    final VariantOutputVcfFormatter formatter = new VariantOutputVcfFormatter(params, "SAMPLE");

    //chr10   2028433 e       A       C:G     11.0    39      4.118   A       3       0.301   C       19      1.900   G       17     1.917
    final VariantSample vs = createSample(Ploidy.DIPLOID, "C:G", false, 11.0 * MathUtils.LOG_10, VariantSample.DeNovoStatus.UNSPECIFIED, null);
    vs.setCoverage(39);
    vs.setCoverageCorrection(4.118);
    vs.setStatisticsString("A\t3\t0.301\tC\t19\t1.900\tG\t17\t1.917");

    //chr10      2028433 .       A       C,G     47      PASS    NS=1;DP=39;AF=0.487,0.436;RE=4.118      GT:GQ:RP:DP:RS  1/2:47:11.0:39:A,3,0.301,C,19,1.900,G,17,1.917
    final VariantLocus locus = new VariantLocus("chr10", 2028432, 2028433, "A", 'N');
    final Variant v = new Variant(locus, vs);
    final String s = formatter.formatCall(v);

    assertEquals("chr10\t2028433\t.\tA\tC,G\t.\tPASS\t.\tGT:DP:RE:GQ:RS\t1/2:39:4.118:110:A,3,0.301,C,19,1.900,G,17,1.917\n", s);
  }

  public void test4() {
    final VariantParams params = new VariantParamsBuilder().maxAmbiguity(0.5).create();
    final VariantOutputVcfFormatter formatter = new VariantOutputVcfFormatter(params, "SAMPLE");

    //chr10   2399391 o       T       i       9.8     35      2.602   T       26      2.601   i       9       0.000
    //something like AT A 1/1

    final VariantSample vs = createSample(Ploidy.DIPLOID, ":", false, 9.8 * MathUtils.LOG_10, VariantSample.DeNovoStatus.UNSPECIFIED, null);
    vs.setCoverage(35);
    vs.setCoverageCorrection(2.602);

    final VariantLocus locus = new VariantLocus("chr10", 2399390, 2399391, "T", 'A');
    final Variant v = new Variant(locus, vs);
    final String s = formatter.formatCall(v);
    assertEquals("chr10\t2399390\t.\tAT\tA\t.\tPASS\t.\tGT:DP:RE:GQ\t1/1:35:2.602:98\n", s);
  }

  public void test5() {
    final VariantParams params = new VariantParamsBuilder().maxAmbiguity(0.5).create();
    final VariantOutputVcfFormatter formatter = new VariantOutputVcfFormatter(params, "SAMPLE");

    //chr10   2402489 e       .t      i:GT    12.8    34      3.428   G       13      1.301   GT      2       0.200   GTT     1    0.100   T       18      1.827
    // X XGT  0/1
    final VariantSample vs = createSample(Ploidy.DIPLOID, ":GT", false, 12.8 * MathUtils.LOG_10, VariantSample.DeNovoStatus.UNSPECIFIED, null);
    vs.setCoverage(34);
    vs.setCoverageCorrection(3.428);

    final VariantLocus locus = new VariantLocus("chr10", 2402489, 2402489, "", 'A');
    final Variant v = new Variant(locus, vs);
    v.setNonIdentityPosterior(Double.POSITIVE_INFINITY);
    final String s = formatter.formatCall(v);
    assertEquals("chr10\t2402489\t.\tA\tAGT\t2147483647.0\tPASS\t.\tGT:DP:RE:GQ\t0/1:34:3.428:128\n", s);
  }

  public void testFormatCall() {
    final VariantParams params = new VariantParamsBuilder().maxAmbiguity(0.5).maxCoverageFilter(new StaticThreshold(20)).create();
    final VariantOutputVcfFormatter formatter = new VariantOutputVcfFormatter(params, "SAMPLE");

    final VariantLocus locus = new VariantLocus("chr10", 82349, 82350, "T", 'A');
    final Variant v = new Variant(locus, createSample(Ploidy.DIPLOID, "A:A", false, 0.0, VariantSample.DeNovoStatus.UNSPECIFIED, null));
    final String s = formatter.formatCall(v);
    assertFalse(s, s.startsWith("#"));
  }

  public void testPloidyCall() {
    final VariantParams params = new VariantParamsBuilder().maxAmbiguity(0.5).maxCoverageFilter(new StaticThreshold(20)).create();
    final VariantOutputVcfFormatter formatter = new VariantOutputVcfFormatter(params, "SAMPLE");
    final VariantLocus locus = new VariantLocus("test", 82349, 82350, "T", 'C');
    final Variant v = new Variant(locus, createSample(Ploidy.HAPLOID, "A", false, 1.0, VariantSample.DeNovoStatus.UNSPECIFIED, null));
    String s = formatter.formatCall(v);
    assertEquals("test\t82350\t.\tT\tA\t.\tPASS\t.\tGT:GQ\t1:6\n", s);
    final Variant v2 = new Variant(locus, createSample(Ploidy.DIPLOID, "A:A", false, 1.0, VariantSample.DeNovoStatus.UNSPECIFIED, null));
    s = formatter.formatCall(v2);
    assertEquals("test\t82350\t.\tT\tA\t.\tPASS\t.\tGT:GQ\t1/1:6\n", s);

    final Variant v3 = new Variant(locus, createSample(Ploidy.HAPLOID, "T", true, null, VariantSample.DeNovoStatus.UNSPECIFIED, null));
    s = formatter.formatCall(v3);
    assertEquals("test\t82350\t.\tT\t.\t.\tPASS\t.\tGT:GQ\t0:3\n", s);
    final Variant v4 = new Variant(locus, createSample(Ploidy.DIPLOID, "T:T", true, null, VariantSample.DeNovoStatus.UNSPECIFIED, null));
    s = formatter.formatCall(v4);
    assertEquals("test\t82350\t.\tT\t.\t.\tPASS\t.\tGT:GQ\t0/0:3\n", s);

    final Variant v5 = new Variant(locus, createSample(Ploidy.HAPLOID, null, true, null, VariantSample.DeNovoStatus.UNSPECIFIED, null));
    s = formatter.formatCall(v5);
    assertEquals("test\t82350\t.\tT\t.\t.\tPASS\t.\tGT\t.\n", s);
    final Variant v6 = new Variant(locus, createSample(Ploidy.DIPLOID, null, true, null, VariantSample.DeNovoStatus.UNSPECIFIED, null));
    s = formatter.formatCall(v6);
    assertEquals("test\t82350\t.\tT\t.\t.\tPASS\t.\tGT\t./.\n", s);
  }

  public void testHeteroComplexInsert() {
    final VariantParams params = new VariantParamsBuilder().maxAmbiguity(0.5).create();
    final VariantOutputVcfFormatter formatter = new VariantOutputVcfFormatter(params, "SAMPLE");

    final VariantLocus locus = new VariantLocus("templateName", 82349, 82349, "", 'C');
    final Variant v = new Variant(locus, createSample(Ploidy.DIPLOID, "C:", false, 18.96, VariantSample.DeNovoStatus.UNSPECIFIED, null));
    v.setIndel(0);
    final String s = formatter.formatCall(v);
    assertEquals("templateName\t82349\t.\tC\tCC\t.\tPASS\t.\tGT:GQ\t1/0:82\n", s);
  }

  public void testComplexX() {
    final VariantParams params = new VariantParamsBuilder().maxAmbiguity(0.5).create();
    final VariantOutputVcfFormatter formatter = new VariantOutputVcfFormatter(params, "SAMPLE");

    final VariantLocus locus = new VariantLocus("templateName", 82349, 82350, "C", 'T');
    final Variant v = new Variant(locus, createSample(Ploidy.DIPLOID, "x", false, 0.0, VariantSample.DeNovoStatus.UNSPECIFIED, null));
    v.addFilter(VariantFilter.FAILED_COMPLEX);
    final String s = formatter.formatCall(v);
    assertEquals("templateName\t82350\t.\tC\t.\t.\tRC\t.\tGT\t./.\n", s);
  }

  public void testNonComplexDelete() {
    final VariantParams params = new VariantParamsBuilder().maxAmbiguity(0.5).create();
    final VariantOutputVcfFormatter formatter = new VariantOutputVcfFormatter(params, "SAMPLE");

    final VariantLocus locus = new VariantLocus("templateName", 82349, 82350, "C", 'C');
    final Variant v = new Variant(locus, createSample(Ploidy.DIPLOID, "T:", false, 0.14, VariantSample.DeNovoStatus.UNSPECIFIED, null));
    final String s = formatter.formatCall(v);
    assertEquals("templateName\t82349\t.\tCC\tCT,C\t.\tPASS\t.\tGT:GQ\t1/2:3\n", s);
  }

  public void testComplexInsert() {
    final VariantParams params = new VariantParamsBuilder().maxAmbiguity(0.5).create();
    final VariantOutputVcfFormatter formatter = new VariantOutputVcfFormatter(params, "SAMPLE");

    final VariantLocus locus = new VariantLocus("templateName", 82349, 82351, "CT", 'C');
    final Variant v = new Variant(locus, createSample(Ploidy.DIPLOID, "CT:AATAA", false, 0.14, VariantSample.DeNovoStatus.UNSPECIFIED, null));
    final String s = formatter.formatCall(v);
    assertEquals("templateName\t82350\t.\tCT\tAATAA\t.\tPASS\t.\tGT:GQ\t0/1:3\n", s);
  }

  public void testComplexHomozygousInsert() {
    final VariantParams params = new VariantParamsBuilder().maxAmbiguity(0.5).create();
    final VariantOutputVcfFormatter formatter = new VariantOutputVcfFormatter(params, "SAMPLE");

    final VariantLocus locus = new VariantLocus("templateName", 82349, 82358, "GGTGTGGAA", 'C');
    final Variant v = new Variant(locus, createSample(Ploidy.DIPLOID, "AGATTGAGTG:AGATTGAGTG", false, 0.14, VariantSample.DeNovoStatus.UNSPECIFIED, null));
    final String s = formatter.formatCall(v);
    assertEquals("templateName\t82350\t.\tGGTGTGGAA\tAGATTGAGTG\t.\tPASS\t.\tGT:GQ\t1/1:3\n", s);
  }

  public void testComplexDelete() {
    final VariantParams params = new VariantParamsBuilder().maxAmbiguity(0.5).create();
    final VariantOutputVcfFormatter formatter = new VariantOutputVcfFormatter(params, "SAMPLE");

    final VariantLocus locus = new VariantLocus("templateName", 82349, 82352, "CCA", 'G');
    final Variant v = new Variant(locus, createSample(Ploidy.DIPLOID, "TC:GC", false, 0.14, VariantSample.DeNovoStatus.UNSPECIFIED, null));
    final String s = formatter.formatCall(v);
    assertEquals("templateName\t82350\t.\tCCA\tTC,GC\t.\tPASS\t.\tGT:GQ\t1/2:3\n", s);
  }

  public void testComplexDeleteAndInsert() {
    final VariantParams params = new VariantParamsBuilder().maxAmbiguity(0.5).create();
    final VariantOutputVcfFormatter formatter = new VariantOutputVcfFormatter(params, "SAMPLE");

    final VariantLocus locus = new VariantLocus("templateName", 82349, 82352, "CCA", 'G');
    final Variant v = new Variant(locus, createSample(Ploidy.DIPLOID, "TC:GCCA", false, 0.14, VariantSample.DeNovoStatus.UNSPECIFIED, null));
    final String s = formatter.formatCall(v);
    assertEquals("templateName\t82350\t.\tCCA\tTC,GCCA\t.\tPASS\t.\tGT:GQ\t1/2:3\n", s);
  }

  public void testComplexDelete2() {
    final VariantParams params = new VariantParamsBuilder().maxAmbiguity(0.5).create();
    final VariantOutputVcfFormatter formatter = new VariantOutputVcfFormatter(params, "SAMPLE");

    final VariantLocus locus = new VariantLocus("templateName", 82349, 82351, "AG", 'C');
    final Variant v = new Variant(locus, createSample(Ploidy.DIPLOID, "AG:", false, 0.14, VariantSample.DeNovoStatus.UNSPECIFIED, null));
    final String s = formatter.formatCall(v);
    assertEquals("templateName\t82349\t.\tCAG\tC\t.\tPASS\t.\tGT:GQ\t0/1:3\n", s);
  }

  public void testComplexInsert2() {
    final VariantParams params = new VariantParamsBuilder().maxAmbiguity(0.5).create();
    final VariantOutputVcfFormatter formatter = new VariantOutputVcfFormatter(params, "SAMPLE");

    final VariantLocus locus = new VariantLocus("templateName", 82349, 82351, "CT", 'C');
    final Variant v = new Variant(locus, createSample(Ploidy.DIPLOID, "TAC:CT", false, 0.14, VariantSample.DeNovoStatus.UNSPECIFIED, null));
    final String s = formatter.formatCall(v);
    assertEquals("templateName\t82350\t.\tCT\tTAC\t.\tPASS\t.\tGT:GQ\t1/0:3\n", s);
  }

  public void testComplexDelete3() {
    final VariantParams params = new VariantParamsBuilder().maxAmbiguity(0.5).create();
    final VariantOutputVcfFormatter formatter = new VariantOutputVcfFormatter(params, "SAMPLE");

    final VariantLocus locus = new VariantLocus("templateName", 82349, 82358, "GGAACGTAA", 'A');
    final Variant v = new Variant(locus, createSample(Ploidy.DIPLOID, "CTTAAT:CTTAAT", false, 0.14, VariantSample.DeNovoStatus.UNSPECIFIED, null));
    final String s = formatter.formatCall(v);
    assertEquals("templateName\t82350\t.\tGGAACGTAA\tCTTAAT\t.\tPASS\t.\tGT:GQ\t1/1:3\n", s);
  }

  public void testComplexCall0Length() {
    final VariantParams params = new VariantParamsBuilder().maxAmbiguity(0.5).create();
    final VariantOutputVcfFormatter formatter = new VariantOutputVcfFormatter(params, "SAMPLE");

    final VariantLocus locus = new VariantLocus("templateName", 82349, 82349, "", 'A');
    final Variant v = new Variant(locus, createSample(Ploidy.DIPLOID, "x", false, 0.0, VariantSample.DeNovoStatus.UNSPECIFIED, null));
    v.setComplexScored();
    v.addFilter(VariantFilter.FAILED_COMPLEX);
    final String s = formatter.formatCall(v);
    assertEquals("templateName\t82349\t.\tA\t.\t.\tRC\tXRX\tGT\t./.\n", s);
  }

  public void testComplexCallNoStats() {
    final VariantParams params = new VariantParamsBuilder().maxAmbiguity(0.5).create();
    final VariantOutputVcfFormatter formatter = new VariantOutputVcfFormatter(params, "SAMPLE");

    final VariantLocus locus = new VariantLocus("templateName", 82349, 82351, "TT", 'A');
    final VariantSample sample = createSample(Ploidy.DIPLOID, "CT:CT", false, 0.0, VariantSample.DeNovoStatus.UNSPECIFIED, null);
    sample.setStatisticsString("blahblah");
    final Variant v = new Variant(locus, sample);
    v.setComplexScored();
    final String s = formatter.formatCall(v);
    assertEquals("templateName\t82350\t.\tTT\tCT\t.\tPASS\tXRX\tGT:GQ\t1/1:3\n", s);
  }

  public void testFilters() {
    final VariantParams params = new VariantParamsBuilder().maxAmbiguity(0.5).maxCoverageFilter(new StaticThreshold(20)).create();
    final VariantOutputVcfFormatter formatter = new VariantOutputVcfFormatter(params, "TEST");

    final VariantLocus locus = new VariantLocus("templateName", 82349, 82350, "T", 'A');
    final Variant v = new Variant(locus, createSample(Ploidy.DIPLOID, "A:A", false, 0.0, VariantSample.DeNovoStatus.UNSPECIFIED, null));
    v.addFilter(VariantFilter.HYPER_COMPLEX);
    String s = formatter.formatCall(v);
    assertEquals("templateName\t82350\t.\tT\tA\t.\tRX\t.\tGT:GQ\t1/1:3\n", s);

    v.addFilter(VariantFilter.COVERAGE);
    s = formatter.formatCall(v);
    assertEquals("templateName\t82350\t.\tT\tA\t.\tOC;RX\tCT=20\tGT:GQ\t1/1:3\n", s);

    v.addFilter(VariantFilter.AMBIGUITY);
    s = formatter.formatCall(v);
    assertEquals("templateName\t82350\t.\tT\tA\t.\tOC;a50.0;RX\tCT=20\tGT:GQ\t1/1:3\n", s);
  }

  public void testErrors() throws IOException {
    final GenomeRelationships grf = RelationshipsFileParser.load(new BufferedReader(new StringReader("original-derived normal cancer contamination=0.30")));
    final VariantParams params = new VariantParamsBuilder().genomeRelationships(grf).maxAmbiguity(0.5).vcfRp(true).create();
    try {
      new VariantOutputVcfFormatter(params);
      fail();
    } catch (final IllegalArgumentException e) {
      assertEquals("No sample names", e.getMessage());
    }
  }

  public void testMultigenomeOutput() throws IOException {
    CommandLine.setCommandArgs("yo", "--ho", "ho");

    final GenomeRelationships grf = RelationshipsFileParser.load(new BufferedReader(new StringReader("original-derived normal cancer contamination=0.30")));
    try (MemoryPrintStream mps = new MemoryPrintStream()) {
      final VariantParams params = new VariantParamsBuilder().genomeRelationships(grf).maxAmbiguity(0.5).vcfRp(true).create();
      final VariantOutputVcfFormatter formatter = new VariantOutputVcfFormatter(params, "normal", "cancer");
      final SAMFileHeader samHeader = new SAMFileHeader();
      SAMReadGroupRecord readGroup = new SAMReadGroupRecord("n");
      readGroup.setPlatform("Illumina");
      readGroup.setSample("normal");
      samHeader.addReadGroup(readGroup);
      readGroup = new SAMReadGroupRecord("c");
      readGroup.setPlatform("Illumina");
      readGroup.setSample("cancer");
      samHeader.addReadGroup(readGroup);
      samHeader.setSequenceDictionary(new SAMSequenceDictionary());
      final SAMSequenceRecord sequenceRecord = new SAMSequenceRecord("chr10", 100000);
      sequenceRecord.setAttribute(SAMSequenceRecord.ASSEMBLY_TAG, "some_assembly_needed");
      sequenceRecord.setAttribute(SAMSequenceRecord.MD5_TAG, "some_hash_string");
      sequenceRecord.setAttribute(SAMSequenceRecord.SPECIES_TAG, "timelord");
      samHeader.getSequenceDictionary().addSequence(sequenceRecord);

      formatter.addExtraInfoFields(EnumSet.of(VcfInfoField.SOMATIC, VcfInfoField.RSS));

      formatter.writeHeader(mps.outputStream(), params, samHeader);

      final String header = mps.toString().replaceAll("##source=.*", "##source=").replaceAll("##RUN-ID=.*", "##RUN-ID=").replaceAll("##fileDate=20.*", "##fileDate=20");
      mNano.check("vovf-multi-header1.vcf", header);

      mps.reset();

      final VariantSample vsNormal = createSample(Ploidy.DIPLOID, "C:C", true, 0.8 * MathUtils.LOG_10, VariantSample.DeNovoStatus.UNSPECIFIED, null, 2.0);
      vsNormal.setCoverage(4);
      vsNormal.setCoverageCorrection(0.4);

      VariantSample vsCancer = createSample(Ploidy.DIPLOID, "T:T", false, 0.9 * MathUtils.LOG_10, VariantSample.DeNovoStatus.UNSPECIFIED, null, 3.0);
      vsCancer.setCoverage(5);
      vsCancer.setCoverageCorrection(0.4);

      final VariantLocus locus = new VariantLocus("chr10", 82349, 82350, "C", 'N');
      final Variant v = new Variant(locus, vsNormal, vsCancer);
      String s = formatter.formatCall(v);
      assertEquals("chr10\t82350\t.\tC\tT\t.\tPASS\tSOMATIC=*;DP=9\tGT:DP:RE:GQ:RP\t0/0:4:0.400:9:0.8\t1/1:5:0.400:10:0.9\n", s);

      vsCancer = createSample(Ploidy.DIPLOID, "C:T", false, 0.9 * MathUtils.LOG_10, VariantSample.DeNovoStatus.UNSPECIFIED, null, 3.0);
      vsCancer.setCoverage(5);
      vsCancer.setCoverageCorrection(0.4);

      final Variant v2 = new Variant(locus, vsNormal, vsCancer);
      v2.setNonIdentityPosterior(1.0);
      v2.setPossibleCause("C:T");
      v2.setPossibleCauseScore(4.0);
      s = formatter.formatCall(v2);
      assertEquals("chr10\t82350\t.\tC\tT\t5.7\tPASS\tSOMATIC=C:T;RSS=1.7;DP=9\tGT:DP:RE:GQ:RP\t0/0:4:0.400:9:0.8\t0/1:5:0.400:10:0.9\n", s);

      vsCancer = createSample(Ploidy.DIPLOID, "C:", false, 0.9 * MathUtils.LOG_10, VariantSample.DeNovoStatus.UNSPECIFIED, null, 3.0);
      vsCancer.setCoverage(5);
      vsCancer.setCoverageCorrection(0.4);

      final Variant v3 = new Variant(locus, vsNormal, vsCancer);
      v3.setNonIdentityPosterior(1.0);
      v3.setPossibleCause("C:");
      v3.setPossibleCauseScore(4.0);
      s = formatter.formatCall(v3);
      assertEquals("chr10\t82349\t.\tNC\tN\t5.7\tPASS\tSOMATIC=NC:N;RSS=1.7;DP=9\tGT:DP:RE:GQ:RP\t0/0:4:0.400:9:0.8\t0/1:5:0.400:10:0.9\n", s);

      vsCancer = createSample(Ploidy.DIPLOID, "C", true, 0.9 * MathUtils.LOG_10, VariantSample.DeNovoStatus.UNSPECIFIED, null, 3.0);
      vsCancer.setCoverage(5);
      vsCancer.setCoverageCorrection(0.4);

      final Variant v4 = new Variant(locus, true, vsNormal, vsCancer);
      v4.setNonIdentityPosterior(1.0);
      v4.setPossibleCause("C:");
      v4.setPossibleCauseScore(4.0);
      s = formatter.formatCall(v4);
      assertEquals("chr10\t82350\t.\tC\t.\t1.4\tPASS\tRSS=1.7;DP=9\tGT:DP:RE:GQ:RP\t0/0:4:0.400:9:0.8\t0/0:5:0.400:10:0.9\n", s);
    }
  }

  public void testDisagree() throws IOException {
    final GenomeRelationships grf = RelationshipsFileParser.load(new BufferedReader(new StringReader(""
        + "parent-child FATHER CHILD\n"
        + "parent-child MOTHER CHILD\n"
        )));

    final VariantParams params = new VariantParamsBuilder().genomeRelationships(grf).maxAmbiguity(0.5).vcfRp(true).create();
    final VariantOutputVcfFormatter formatter = new VariantOutputVcfFormatter(params, "FATHER", "MOTHER", "CHILD");

    //disagree
    final VariantLocus locus = new VariantLocus("chrY", 82349, 82350, "C", 'N');
    final VariantSample vs = createSample(Ploidy.HAPLOID, "C", false, null, VariantSample.DeNovoStatus.UNSPECIFIED, null, 0.05);
    vs.setAmbiguityRatio(3.0);
    vs.setCoverage(5, 0.2587);

    //no disagree
    final VariantSample vs2 = createSample(Ploidy.HAPLOID, "C", false, null, VariantSample.DeNovoStatus.UNSPECIFIED, null, 0.05);
    vs2.setAmbiguityRatio(1.0);
    vs2.setCoverage(5, 0.2587);

    final Variant v = new Variant(locus, vs, vs2, null);
    final String actual = formatter.formatCall(v);
    assertEquals("chrY" + TAB + "82350" + TAB + "." + TAB + "C" + TAB + "." + TAB + "." + TAB + "PASS" + TAB + "DP=10" + TAB + "GT:DP:RE:AR:GQ:RP" + TAB + "0:5:0.259:3.000:3:0.0" + TAB + "0:5:0.259:1.000:3:0.0" + TAB + ".\n", actual);
  }


  public void testFamilyOutput() throws IOException {
    final GenomeRelationships grf = RelationshipsFileParser.load(new BufferedReader(new StringReader(""
        + "parent-child FATHER CHILD\n"
        + "parent-child MOTHER CHILD\n"
        )));

    final VariantParams params = new VariantParamsBuilder().genomeRelationships(grf).maxAmbiguity(0.5).vcfRp(true).create();
    final VariantOutputVcfFormatter formatter = new VariantOutputVcfFormatter(params, "FATHER", "MOTHER", "CHILD");
    final SAMFileHeader samHeader = new SAMFileHeader();
    SAMReadGroupRecord readGroup = new SAMReadGroupRecord("F");
    readGroup.setPlatform("Illumina");
    readGroup.setSample("FATHER");
    samHeader.addReadGroup(readGroup);
    readGroup = new SAMReadGroupRecord("M");
    readGroup.setPlatform("Illumina");
    readGroup.setSample("MOTHER");
    samHeader.addReadGroup(readGroup);
    readGroup = new SAMReadGroupRecord("C");
    readGroup.setPlatform("Illumina");
    readGroup.setSample("CHILD");
    samHeader.addReadGroup(readGroup);
    samHeader.setSequenceDictionary(new SAMSequenceDictionary());
    final SAMSequenceRecord sequenceRecord = new SAMSequenceRecord("chr10", 100000);
    sequenceRecord.setAttribute(SAMSequenceRecord.ASSEMBLY_TAG, "some_assembly_needed");
    sequenceRecord.setAttribute(SAMSequenceRecord.MD5_TAG, "some_hash_string");
    sequenceRecord.setAttribute(SAMSequenceRecord.SPECIES_TAG, "timelord");
    samHeader.getSequenceDictionary().addSequence(sequenceRecord);

    //haploid
    final VariantLocus locus = new VariantLocus("chrY", 82349, 82350, "C", 'N');
    VariantSample vs = createSample(Ploidy.HAPLOID, "C", false, null, VariantSample.DeNovoStatus.UNSPECIFIED, null, 0.05);
    vs.setAmbiguityRatio(3.0);
    vs.setCoverage(5, 0.2587);
    Variant v = new Variant(locus, vs, null, null);
    String actual = formatter.formatCall(v);
    assertEquals("chrY" + TAB + "82350" + TAB + "." + TAB + "C" + TAB + "." + TAB + "." + TAB + "PASS" + TAB + "DP=5" + TAB + "GT:DP:RE:AR:GQ:RP" + TAB + "0:5:0.259:3.000:3:0.0" + TAB + "." + TAB + ".\n", actual);

    //diploid homozygous
    vs = createSample(Ploidy.DIPLOID, "C:C", false, null, VariantSample.DeNovoStatus.UNSPECIFIED, null, 0.05);
    vs.setAmbiguityRatio(3.0);
    vs.setCoverage(5, 0.2587);
    v = new Variant(locus, vs, null, null);
    actual = formatter.formatCall(v);
    assertEquals("chrY" + TAB + "82350" + TAB + "." + TAB + "C" + TAB + "." + TAB + "." + TAB + "PASS" + TAB + "DP=5" + TAB + "GT:DP:RE:AR:GQ:RP" + TAB + "0/0:5:0.259:3.000:3:0.0" + TAB + "." + TAB + ".\n", actual);

    //diploid het
    vs = createSample(Ploidy.DIPLOID, "C:C", false, null, VariantSample.DeNovoStatus.UNSPECIFIED, null, 0.05);
    vs.setAmbiguityRatio(3.0);
    vs.setCoverage(5, 0.2587);
    v = new Variant(locus, vs, null, null);
    actual = formatter.formatCall(v);
    assertEquals("chrY" + TAB + "82350" + TAB + "." + TAB + "C" + TAB + "." + TAB + "." + TAB + "PASS" + TAB + "DP=5" + TAB + "GT:DP:RE:AR:GQ:RP" + TAB + "0/0:5:0.259:3.000:3:0.0" + TAB + "." + TAB + ".\n", actual);

    //diploid both het indel
    vs = createSample(Ploidy.DIPLOID, "C:", false, null, VariantSample.DeNovoStatus.UNSPECIFIED, null, 0.05);
    vs.setAmbiguityRatio(3.0);
    vs.setCoverage(5, 0.2587);
    v = new Variant(locus, vs, null, null);
    actual = formatter.formatCall(v);
    assertEquals("chrY" + TAB + "82349" + TAB + "." + TAB + "NC" + TAB + "N" + TAB + "." + TAB + "PASS" + TAB + "DP=5" + TAB + "GT:DP:RE:AR:GQ:RP" + TAB + "0/1:5:0.259:3.000:3:0.0" + TAB + "." + TAB + ".\n", actual);

    //diploid het indel
    vs = createSample(Ploidy.DIPLOID, "C:C", false, null, VariantSample.DeNovoStatus.UNSPECIFIED, null, 0.05);
    vs.setAmbiguityRatio(3.0);
    vs.setCoverage(5, 0.2587);
    v = new Variant(locus, vs, null, null);
    actual = formatter.formatCall(v);
    assertEquals("chrY" + TAB + "82350" + TAB + "." + TAB + "C" + TAB + "." + TAB + "." + TAB + "PASS" + TAB + "DP=5" + TAB + "GT:DP:RE:AR:GQ:RP" + TAB + "0/0:5:0.259:3.000:3:0.0" + TAB + "." + TAB + ".\n", actual);
  }

  public void testFamilyOutputDenovo() throws IOException {
    final GenomeRelationships grf = RelationshipsFileParser.load(new BufferedReader(new StringReader(""
        + "parent-child FATHER CHILD\n"
        + "parent-child MOTHER CHILD\n"
        )));

    final VariantParams params = new VariantParamsBuilder().genomeRelationships(grf).maxAmbiguity(0.5).vcfRp(true).create();
    final VariantOutputVcfFormatter formatter = new VariantOutputVcfFormatter(params, "FATHER", "MOTHER", "CHILD");
    final SAMFileHeader samHeader = new SAMFileHeader();
    SAMReadGroupRecord readGroup = new SAMReadGroupRecord("F");
    readGroup.setPlatform("Illumina");
    readGroup.setSample("FATHER");
    samHeader.addReadGroup(readGroup);
    readGroup = new SAMReadGroupRecord("M");
    readGroup.setPlatform("Illumina");
    readGroup.setSample("MOTHER");
    samHeader.addReadGroup(readGroup);
    readGroup = new SAMReadGroupRecord("C");
    readGroup.setPlatform("Illumina");
    readGroup.setSample("CHILD");
    samHeader.addReadGroup(readGroup);
    samHeader.setSequenceDictionary(new SAMSequenceDictionary());
    final SAMSequenceRecord sequenceRecord = new SAMSequenceRecord("chr10", 100000);
    sequenceRecord.setAttribute(SAMSequenceRecord.ASSEMBLY_TAG, "some_assembly_needed");
    sequenceRecord.setAttribute(SAMSequenceRecord.MD5_TAG, "some_hash_string");
    sequenceRecord.setAttribute(SAMSequenceRecord.SPECIES_TAG, "timelord");
    samHeader.getSequenceDictionary().addSequence(sequenceRecord);

    final VariantLocus locus = new VariantLocus("chr1", 82349, 82350, "C", 'A');
    final VariantSample father = createSample(Ploidy.HAPLOID, "C", false, null, VariantSample.DeNovoStatus.NOT_DE_NOVO, 0.0, 0.05);
    father.setAmbiguityRatio(3.0);
    father.setCoverage(5, 0.2587);
    final VariantSample mother = createSample(Ploidy.HAPLOID, "G", false, null, VariantSample.DeNovoStatus.NOT_DE_NOVO, 0.0, 0.05);
    mother.setAmbiguityRatio(3.0);
    mother.setCoverage(5, 0.2587);

    final VariantSample son = createSample(Ploidy.HAPLOID, "G:T", false, null, VariantSample.DeNovoStatus.IS_DE_NOVO, 10.0, 0.05);
    son.setAmbiguityRatio(3.0);
    son.setCoverage(5, 0.2587);

    final Variant v = new Variant(locus, father, mother, son);
    final String actual = formatter.formatCall(v);
    //chr1  82350 . C G,G:T . PASS  DP=15 GT:DP:RE:AR:AB:RQ:GQ:RP:DN
    //0:5:0.259:3.000:0.200:3.1:3:0.0:N
    //1:5:0.259:3.000:0.200:3.1:3:0.0:N
    //2:5:0.259:3.000:0.200:3.1:3:0.0:Y
    assertEquals("chr1" + TAB + "82350" + TAB + "." + TAB + "C" + TAB + "G,G:T" + TAB + "." + TAB + "PASS" + TAB + "DP=15" + TAB + "GT:DP:RE:AR:GQ:RP:DN:DNP" + TAB
        + "0:5:0.259:3.000:3:0.0:N:3" + TAB + "1:5:0.259:3.000:3:0.0:N:3"
        + TAB + "2:5:0.259:3.000:3:0.0:Y:43\n", actual);
  }

  public void testDiseasedFamilyOutput() throws IOException {
    CommandLine.setCommandArgs("yo", "--ho", "ho");
    final GenomeRelationships grf = RelationshipsFileParser.load(new BufferedReader(new StringReader(""
        + "genome FATHER disease=true\n"
        + "genome CHILD disease=true\n"
        + "parent-child FATHER CHILD\n"
        + "parent-child MOTHER CHILD\n"
        )));
    try (MemoryPrintStream mps = new MemoryPrintStream()) {
      final VariantParams params = new VariantParamsBuilder().genomeRelationships(grf).maxAmbiguity(0.5).vcfRp(true).create();
      final VariantOutputVcfFormatter formatter = new VariantOutputVcfFormatter(params, "FATHER", "MOTHER", "CHILD");
      final SAMFileHeader samHeader = new SAMFileHeader();
      SAMReadGroupRecord readGroup = new SAMReadGroupRecord("F");
      readGroup.setPlatform("Illumina");
      readGroup.setSample("FATHER");
      samHeader.addReadGroup(readGroup);
      readGroup = new SAMReadGroupRecord("M");
      readGroup.setPlatform("Illumina");
      readGroup.setSample("MOTHER");
      samHeader.addReadGroup(readGroup);
      readGroup = new SAMReadGroupRecord("C");
      readGroup.setPlatform("Illumina");
      readGroup.setSample("CHILD");
      samHeader.addReadGroup(readGroup);
      samHeader.setSequenceDictionary(new SAMSequenceDictionary());
      final SAMSequenceRecord sequenceRecord = new SAMSequenceRecord("chr10", 100000);
      sequenceRecord.setAttribute(SAMSequenceRecord.ASSEMBLY_TAG, "some_assembly_needed");
      sequenceRecord.setAttribute(SAMSequenceRecord.MD5_TAG, "some_hash_string");
      sequenceRecord.setAttribute(SAMSequenceRecord.SPECIES_TAG, "timelord");
      samHeader.getSequenceDictionary().addSequence(sequenceRecord);
      formatter.addExtraInfoFields(EnumSet.of(VcfInfoField.DISEASE, VcfInfoField.RDS));
      formatter.writeHeader(mps.outputStream(), params, samHeader);

      final String header = mps.toString().replaceAll("##source=.*", "##source=").replaceAll("##RUN-ID=.*", "##RUN-ID=").replaceAll("##fileDate=20.*", "##fileDate=20");
      mNano.check("vovf-disease-header1.vcf", header);
    }
  }

}
