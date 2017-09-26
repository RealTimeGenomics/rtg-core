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

import com.rtg.launcher.AbstractNanoTest;
import com.rtg.reference.Ploidy;
import com.rtg.relation.GenomeRelationships;
import com.rtg.relation.RelationshipsFileParser;
import com.rtg.util.MathUtils;
import com.rtg.util.TestUtils;
import com.rtg.util.cli.CommandLine;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.MemoryPrintStream;
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

/**
 */
public class VariantOutputVcfFormatterTest extends AbstractNanoTest {

  private static final String TAB = "\t";

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

        final String header = mps.toString();
        mNano.check("vovf-single-header1.vcf", TestUtils.sanitizeVcfHeader(header));
      }
    }
  }

  public void testVariants() throws IOException {
    final VariantParams params = new VariantParamsBuilder().maxAmbiguity(0.5).create();
    final VariantOutputVcfFormatter formatter = new VariantOutputVcfFormatter(params, "SAMPLE");

    final Variant[] v = {
      variant1(),
      variant2(),
      variant3(),
      variant4(),
      variant5(),
      variant6(),
    };
    final StringBuilder sb = new StringBuilder();
    for (Variant v1 : v) {
      sb.append(formatter.formatCall(v1));
    }
    mNano.check("vovf-test-several.vcf", sb.toString());
  }

  protected Variant variant1() {
    final VariantSample vs = createSample(Ploidy.DIPLOID, "C:T", false, 9.5 * MathUtils.LOG_10, VariantSample.DeNovoStatus.NOT_DE_NOVO, 0.0);
    vs.setCoverage(38);
    vs.setCoverageCorrection(4.622);
    return new Variant(new VariantLocus("chr10", 112108, 112109, "C", 'N'), vs);
  }

  protected Variant variant2() {
    final VariantSample vs = createSample(Ploidy.DIPLOID, "C:G", false, 16.4 * MathUtils.LOG_10, VariantSample.DeNovoStatus.UNSPECIFIED, null);
    vs.setCoverage(81);
    vs.setCoverageCorrection(15.368);
    return new Variant(new VariantLocus("chr10", 115007, 115008, "G", 'N'), vs);
  }

  protected Variant variant3() {
    final VariantSample vs = createSample(Ploidy.DIPLOID, "C:G", false, 11.0 * MathUtils.LOG_10, VariantSample.DeNovoStatus.UNSPECIFIED, null);
    vs.setCoverage(39);
    vs.setCoverageCorrection(4.118);
    vs.setStatisticsString("A\t3\t0.301\tC\t19\t1.900\tG\t17\t1.917");
    return new Variant(new VariantLocus("chr10", 2028432, 2028433, "A", 'N'), vs);
  }

  protected Variant variant4() {
    final VariantSample vs = createSample(Ploidy.DIPLOID, ":", false, 9.8 * MathUtils.LOG_10, VariantSample.DeNovoStatus.UNSPECIFIED, null);
    vs.setCoverage(35);
    vs.setCoverageCorrection(2.602);
    return new Variant(new VariantLocus("chr10", 2399390, 2399391, "T", 'A'), vs);
  }

  protected Variant variant5() {
    final VariantSample vs = createSample(Ploidy.DIPLOID, ":GT", false, 12.8 * MathUtils.LOG_10, VariantSample.DeNovoStatus.UNSPECIFIED, null);
    vs.setCoverage(34);
    vs.setCoverageCorrection(3.428);
    final Variant v = new Variant(new VariantLocus("chr10", 2402489, 2402489, "", 'A'), vs);
    v.setNonIdentityPosterior(Double.POSITIVE_INFINITY);
    return v;
  }

  protected Variant variant6() {
    final VariantSample vs = createSample(Ploidy.DIPLOID, "T:T", false, 0.8 * MathUtils.LOG_10, VariantSample.DeNovoStatus.NOT_DE_NOVO, 0.0);
    vs.setCoverage(4);
    vs.setCoverageCorrection(0.4);
    final Variant v = new Variant(new VariantLocus("chr10", 82349, 82350, "C", 'N'), vs);
    v.setNonIdentityPosterior(2.0);
    return v;
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

      formatter.addExtraFormatFields(EnumSet.of(VcfFormatField.SSC, VcfFormatField.SS));

      formatter.writeHeader(mps.outputStream(), params, samHeader);

      final String header = mps.toString();
      TestUtils.containsAll(header, "yo --ho ho");
      mNano.check("vovf-multi-header1.vcf", TestUtils.sanitizeVcfHeader(header));

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
      assertEquals("chr10\t82350\t.\tC\tT\t.\tPASS\tDP=9\tGT:DP:RE:GQ:RP\t0/0:4:0.400:9:0.8\t1/1:5:0.400:10:0.9\n", s);

      vsCancer = createSample(Ploidy.DIPLOID, "C:T", false, 0.9 * MathUtils.LOG_10, VariantSample.DeNovoStatus.IS_DE_NOVO, 4.0, null);
      vsCancer.setCoverage(5);
      vsCancer.setCoverageCorrection(0.4);

      final Variant v2 = new Variant(locus, vsNormal, vsCancer);
      v2.setNonIdentityPosterior(1.0);
      v2.setPossibleCause("C:T");
      v2.setNormalCancerScore(4.0);
      s = formatter.formatCall(v2);
      assertEquals("chr10\t82350\t.\tC\tT\t5.7\tPASS\tDP=9\tGT:DP:RE:GQ:RP:SSC:SS\t0/0:4:0.400:9:0.8\t0/1:5:0.400:10:0.9:1.7:2\n", s);

      vsCancer = createSample(Ploidy.DIPLOID, "C:", false, 0.9 * MathUtils.LOG_10, VariantSample.DeNovoStatus.IS_DE_NOVO, 4.0, null);
      vsCancer.setCoverage(5);
      vsCancer.setCoverageCorrection(0.4);

      final Variant v3 = new Variant(locus, vsNormal, vsCancer);
      v3.setNonIdentityPosterior(1.0);
      v3.setPossibleCause("C:");
      v3.setNormalCancerScore(4.0);
      s = formatter.formatCall(v3);
      assertEquals("chr10\t82349\t.\tNC\tN\t5.7\tPASS\tDP=9\tGT:DP:RE:GQ:RP:SSC:SS\t0/0:4:0.400:9:0.8\t0/1:5:0.400:10:0.9:1.7:2\n", s);

      vsCancer = createSample(Ploidy.DIPLOID, "C", true, 0.9 * MathUtils.LOG_10, VariantSample.DeNovoStatus.NOT_DE_NOVO, 4.0, null);
      vsCancer.setCoverage(5);
      vsCancer.setCoverageCorrection(0.4);

      final Variant v4 = new Variant(locus, vsNormal, vsCancer);
      v4.setNonIdentityPosterior(1.0);
      v4.setPossibleCause("C:");
      s = formatter.formatCall(v4);
      assertEquals("chr10\t82350\t.\tC\t.\t1.4\tPASS\tDP=9\tGT:DP:RE:GQ:RP:SSC:SS\t0/0:4:0.400:9:0.8\t0/0:5:0.400:10:0.9:1.7:0\n", s);
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
    formatter.addExtraFormatFields(EnumSet.of(VcfFormatField.RQ, VcfFormatField.DN, VcfFormatField.DNP));
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

    final VariantSample son = createSample(Ploidy.HAPLOID, "T", false, null, VariantSample.DeNovoStatus.IS_DE_NOVO, 10.0, 0.05);
    son.setAmbiguityRatio(3.0);
    son.setCoverage(5, 0.2587);

    final Variant v = new Variant(locus, father, mother, son);
    final String actual = formatter.formatCall(v);
    //chr1  82350 . C G,G:T . PASS  DP=15 GT:DP:RE:AR:AB:RQ:GQ:RP:DN
    //0:5:0.259:3.000:0.200:3.1:3:0.0:N
    //1:5:0.259:3.000:0.200:3.1:3:0.0:N
    //2:5:0.259:3.000:0.200:3.1:3:0.0:Y
    assertEquals("chr1" + TAB + "82350" + TAB + "." + TAB + "C" + TAB + "G,T" + TAB + "." + TAB + "PASS" + TAB + "DP=15" + TAB + "GT:DP:RE:AR:GQ:RP:DN:DNP" + TAB
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
      formatter.addExtraFormatFields(EnumSet.of(VcfFormatField.RQ, VcfFormatField.DN, VcfFormatField.DNP));
      formatter.writeHeader(mps.outputStream(), params, samHeader);

      final String header = mps.toString();
      mNano.check("vovf-disease-header1.vcf", TestUtils.sanitizeVcfHeader(header));
    }
  }

}
