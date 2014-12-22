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

package com.rtg.variant.eval;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.StringReader;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Map.Entry;
import java.util.TreeMap;

import com.rtg.launcher.OutputParams;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.reader.SdfId;
import com.rtg.tabix.TabixIndexer;
import com.rtg.tabix.UnindexableDataException;
import com.rtg.util.MathUtils;
import com.rtg.util.Pair;
import com.rtg.util.StringUtils;
import com.rtg.util.TestUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.util.test.FileHelper;
import com.rtg.util.test.NanoRegression;
import com.rtg.variant.PosteriorUtils;
import com.rtg.variant.eval.SequenceEvaluator.VariantPositionComparator;
import com.rtg.vcf.VcfReader;
import com.rtg.vcf.VcfRecord;
import com.rtg.vcf.VcfUtils;
import com.rtg.vcf.header.VcfHeader;

import junit.framework.TestCase;

/**
 */
public class VcfEvalTaskTest extends TestCase {
  private File mDir = null;
  NanoRegression mNano = null;
  @Override
  public void setUp() throws IOException {
    mDir = FileUtils.createTempDir("MutationEval", "mDir");
    Diagnostic.setLogStream(TestUtils.getNullPrintStream());
    mNano = new NanoRegression(this.getClass());
  }

  @Override
  public void tearDown() throws IOException {
    FileHelper.deleteAll(mDir);
    mDir = null;
    Diagnostic.setLogStream();
    try {
      mNano.finish();
    } finally {
      mNano = null;
    }
  }

  private static final String TEMPLATE = ">seq" + StringUtils.LS
    + "ACAGTCACGGTACGTACGTACGTACGT" + StringUtils.LS;

  private static final String[] TP = {"seq 3 . A G 0.0 PASS . GT 1/1"};
  private static final String[] TP_OUT = TP;
  private static final String[] FP = {"seq 6 . C T 0.0 PASS . GT 0/1"};
  private static final String[] FP_OUT =  FP;
  private static final String[] FN = {"seq 9 . G T 0.0 PASS . GT 1/1"};
  private static final String[] FN_OUT = FN;

  /** vcf header */
  private static final String CALLS_HEADER = VcfHeader.MINIMAL_HEADER + "\tRTG";
  private static final String MUTATIONS_HEADER = VcfHeader.MINIMAL_HEADER + "\tRTG";

  public void test() throws IOException, UnindexableDataException {
    check(TP, FP, FN, TP_OUT, FP_OUT, FN_OUT, VcfUtils.FORMAT_GENOTYPE_QUALITY);
  }

  private void check(String[] tp, String[] fp, String[] fn, String[] tpOut, String[] fpOut, String[] fnOut, String sortField) throws IOException, UnindexableDataException {
    check(TEMPLATE, tp, fp, fn, tpOut, fpOut, fnOut, true, sortField);
    check(TEMPLATE, tp, fp, fn, tpOut, fpOut, fnOut, false, sortField);
  }
  private void check(String ref, String[] tp, String[] fp, String[] fn, String[] tpOut, String[] fpOut, String[] fnOut, boolean zip, String sortField) throws IOException, UnindexableDataException {
    createInput(tp, fp, fn);

    final File calls = new File(mDir, "calls.vcf.gz");
    final File mutations = new File(mDir, "mutations.vcf.gz");
    final File out = FileUtils.createTempDir("out", zip ? "zip" : "notZipped", mDir);
    final File template = new File(mDir, "template");
    //System.err.println("baseline \n" + FileUtils.fileToString(mutations));
    //System.err.println("calls \n" + FileUtils.fileToString(calls));
    ReaderTestUtils.getReaderDNA(ref, template, null).close();
    final VcfEvalParams params = VcfEvalParams.builder().baseLineFile(mutations).callsFile(calls)
        .templateFile(template).outputParams(new OutputParams(out, false, zip)).scoreField(sortField).create();
    VcfEvalTask.evaluateCalls(params);
    final String zipSuffix = zip ? ".gz" : "";
    final File tpOutFile = new File(out, "tp.vcf" + zipSuffix);
    final File fpOutFile = new File(out, "fp.vcf" + zipSuffix);
    final File fnOutFile = new File(out, "fn.vcf" + zipSuffix);
    final String tpResults = fileToString(tpOutFile, zip);
    final String fpResults = fileToString(fpOutFile, zip);
    final String fnResults = fileToString(fnOutFile, zip);
    final String header = "##fileformat=VCFv4.1\n" + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tRTG";
    if (zip) {
      assertTrue(new File(out, "tp.vcf" + zipSuffix + TabixIndexer.TABIX_EXTENSION).exists());
      assertTrue(new File(out, "fp.vcf" + zipSuffix + TabixIndexer.TABIX_EXTENSION).exists());
      assertTrue(new File(out, "fn.vcf" + zipSuffix + TabixIndexer.TABIX_EXTENSION).exists());
    }

    //System.err.println("TP \n" + tpResults);
    //System.err.println("FP \n" + fpResults);
    //System.err.println("FN \n" + fnResults);
    TestUtils.containsAll(tpResults, tpOut);
    TestUtils.containsAll(tpResults, header.replaceAll("\t", " "));
    TestUtils.containsAll(fpResults, fpOut);
    TestUtils.containsAll(fpResults, header.replaceAll("\t", " "));
    TestUtils.containsAll(fnResults, fnOut);
    TestUtils.containsAll(fnResults, header.replaceAll("\t", " "));
    final File phase = new File(out, "phasing.txt");
    assertTrue(phase.exists());
    final String phasing = FileUtils.fileToString(phase);
    TestUtils.containsAll(phasing, "Correct phasings: ", "Incorrect phasings: ", "Unresolvable phasings: ");
  }
  private String fileToString(File f, boolean zip) throws IOException {
    final String result;
    if (zip) {
      result = FileHelper.gzFileToString(f);
    } else {
      result = FileUtils.fileToString(f);
    }
    return result.replaceAll("\t", " ");
  }

  private void createInput(String[] tp, String[] fp, String[] fn) throws IOException, UnindexableDataException {
    final File calls = new File(mDir, "calls.vcf.gz");
    final File mutations = new File(mDir, "mutations.vcf.gz");
    final TreeMap<DetectedVariant, String> callList = new TreeMap<>();
    final TreeMap<DetectedVariant, String> mutationList = new TreeMap<>();
    for (final String var : tp) {
      final VcfRecord rec = VcfReader.vcfLineToRecord(var.replaceAll(" ", "\t"));
      callList.put(new DetectedVariant(rec, 0, RocSortValueExtractor.NULL_EXTRACTOR, false), rec.toString());
      mutationList.put(new DetectedVariant(rec, 0, RocSortValueExtractor.NULL_EXTRACTOR, false), rec.toString());
    }
    for (final String var : fp) {
      final VcfRecord rec = VcfReader.vcfLineToRecord(var.replaceAll(" ", "\t"));
      callList.put(new DetectedVariant(rec, 0, RocSortValueExtractor.NULL_EXTRACTOR, false), rec.toString());
    }
    for (final String var : fn) {
      final VcfRecord rec = VcfReader.vcfLineToRecord(var.replaceAll(" ", "\t"));
      mutationList.put(new DetectedVariant(rec, 0, RocSortValueExtractor.NULL_EXTRACTOR, false), rec.toString());
    }
    try (BufferedWriter callOut = new BufferedWriter(new OutputStreamWriter(FileUtils.createOutputStream(calls, true, false)))) {
      callOut.write(CALLS_HEADER.replaceAll(" ", "\t") + StringUtils.LS);
      for (final Entry<DetectedVariant, String> var : callList.entrySet()) {
        callOut.write(var.getValue() + "\n");
      }
    }
    new TabixIndexer(calls).saveVcfIndex();
    try (BufferedWriter mutOut = new BufferedWriter(new OutputStreamWriter(FileUtils.createOutputStream(mutations, true, false)))) {
      mutOut.write(MUTATIONS_HEADER.replaceAll(" ", "\t") + StringUtils.LS);
      for (final Entry<DetectedVariant, String> var : mutationList.entrySet()) {
        mutOut.write(var.getValue() + "\n");
      }
    }
    new TabixIndexer(mutations).saveVcfIndex();

  }
  //   12 3456 789 012 34567890 123 456789
  //   01 234 5 67 890 12345 6789 012 3456
  //  "AC AGT C AC GGT ACGTA CGTA CGT ACGT"
  //   AC GGT C AC TGT
  //   AC GGT T AC GGT

  private static final String[] TP_LARGE = {
      "seq 3 . A G 5.0 PASS . GT:GQ 1/1:" + PosteriorUtils.phredIfy(5.0 * MathUtils.LOG_10)
    , "seq 12 . A C 3.0 PASS . GT:GQ 0/1:" + PosteriorUtils.phredIfy(3.0 * MathUtils.LOG_10)
    , "seq 14 . G GA,T 1.0 PASS . GT:GQ 1/2:" + PosteriorUtils.phredIfy(1.0 * MathUtils.LOG_10)
  };
  private static final String[] FP_LARGE = {
      "seq 6 . C T 4.0 PASS . GT:GQ 1/0:" + PosteriorUtils.phredIfy(4.0 * MathUtils.LOG_10)
    , "seq 23 . T G 0.0 PASS . GT 1/1"
    , "seq 26 . G C 2.0 PASS . GT:GQ 0/1:" + PosteriorUtils.phredIfy(2.0 * MathUtils.LOG_10)
  };
  private static final String[] FN_LARGE = {
       "seq 9 . G T 0.0 PASS . GT 1/1"
    , "seq 16 . A T 0.0 PASS . GT 0/1"
    , "seq 20 . A T 0.0 PASS . GT 1/1"
  };

  public void testLarge() throws IOException, UnindexableDataException {
    check(TP_LARGE, FP_LARGE, FN_LARGE, TP_LARGE, FP_LARGE, FN_LARGE, VcfUtils.FORMAT_GENOTYPE_QUALITY);
  }


  public void testROC() throws IOException, UnindexableDataException {
    checkRoc("testroc", TEMPLATE, TP_LARGE, FP_LARGE, FN_LARGE);
  }

  private void checkRoc(String label, String template, String[] tp, String[] fp, String[] fn) throws IOException, UnindexableDataException {
    checkRoc(label, template, tp, fp, fn, true);
  }

  private void checkRoc(String label, String template, String[] tp, String[] fp, String[] fn, boolean checktotal) throws IOException, UnindexableDataException {
    checkRoc(label, template, tp, fp, fn, checktotal, true);
    checkRoc(label, template, tp, fp, fn, checktotal, false);
  }
  private void checkRoc(String label, String template, String[] tp, String[] fp, String[] fn, boolean checktotal, boolean rtgStats) throws IOException, UnindexableDataException {
    createInput(tp, fp, fn);

    final File calls = new File(mDir, "calls.vcf.gz");
    final File mutations = new File(mDir, "mutations.vcf.gz");
    final File out = FileUtils.createTempDir("out", "", mDir);
    final File genome = new File(mDir, "template");
    ReaderTestUtils.getReaderDNA(template, genome, null).close();
    final VcfEvalParams params = VcfEvalParams.builder().baseLineFile(mutations).callsFile(calls)
        .templateFile(genome).outputParams(new OutputParams(out, false, false))
        .rtgStats(rtgStats)
        .create();
    VcfEvalTask.evaluateCalls(params);
    final int tpCount = tp.length;
    final int fnCount = fn.length;

    checkRocResults(label + "-weighted.tsv", new File(out, VcfEvalTask.FULL_ROC_FILE), checktotal, tpCount, fnCount);
    checkRocResults(label + "-homo.tsv", new File(out, VcfEvalTask.HOMOZYGOUS_FILE), checktotal, tpCount, fnCount);
    checkRocResults(label + "-hetero.tsv", new File(out, VcfEvalTask.HETEROZYGOUS_FILE), checktotal, tpCount, fnCount);
    if (rtgStats) {
      checkRocResults(label + "-simple.tsv", new File(out, VcfEvalTask.SIMPLE_FILE), checktotal, tpCount, fnCount);
      checkRocResults(label + "-complex.tsv", new File(out, VcfEvalTask.COMPLEX_FILE), checktotal, tpCount, fnCount);
    } else {
      assertFalse(new File(out, VcfEvalTask.SIMPLE_FILE).exists());
      assertFalse(new File(out, VcfEvalTask.COMPLEX_FILE).exists());
    }

    final VcfEvalParams paramsrev = VcfEvalParams.builder().baseLineFile(mutations).callsFile(calls)
        .templateFile(genome).outputParams(new OutputParams(out, false, false))
        .sortOrder(RocSortOrder.ASCENDING)
        .rtgStats(rtgStats)
        .create();
    VcfEvalTask.evaluateCalls(paramsrev);
    checkRocResults(label + "-weighted-rev.tsv", new File(out, VcfEvalTask.FULL_ROC_FILE), checktotal, tpCount, fnCount);
  }

  private void checkRocResults(String label, final File out, boolean checktotal, final int tpCount, final int fnCount) throws IOException {
    final String roc = FileUtils.fileToString(out);
    //System.err.println("ROC\n" + roc);
    final String[] homoLines = roc.split(StringUtils.LS);
    if (checktotal) {
      assertEquals("#total baseline variants: " + (tpCount + fnCount), homoLines[0]);
    } else {
      assertTrue(homoLines[0].startsWith("#total baseline variants: "));
    }
    assertTrue(homoLines[1].startsWith("#score\t"));
    mNano.check(label, roc);
  }

  private static final String[] EMPTY = {};
  private static final String[] FP_TRICKY = {
      "seq 1 . A T 1 PASS . GT:GQ 1/1:" + PosteriorUtils.phredIfy(1 * MathUtils.LOG_10)
    , "seq 2 . A AT 9 PASS . GT:GQ 1/1:" + PosteriorUtils.phredIfy(9 * MathUtils.LOG_10)
    , "seq 5 . GT G 2 PASS . GT:GQ 1/1:" + PosteriorUtils.phredIfy(2 * MathUtils.LOG_10)
    , "seq 10 . CGT AGA 3 PASS . GT:GQ 1/1:" + PosteriorUtils.phredIfy(3 * MathUtils.LOG_10)
    , "seq 14 . C A 8 PASS . GT:GQ 1/1:" + PosteriorUtils.phredIfy(8 * MathUtils.LOG_10)
    , "seq 16 . T A 4 PASS . GT:GQ 1/1:" + PosteriorUtils.phredIfy(4 * MathUtils.LOG_10)
    , "seq 20 . C A 10 PASS . GT:GQ 0/1:" + PosteriorUtils.phredIfy(10 * MathUtils.LOG_10)
  };

  private static final String TRICKY_TEMPLATE = ">seq" + StringUtils.LS
      + "ACTTTCCCACGTACGTCCTCT" + StringUtils.LS;

  private static final String[] FN_TRICKY = {
      "seq 5 . TC T 0.0 PASS . GT 1/1"
    , "seq 7 . C CC 0.0 PASS . GT 1/1"
    , "seq 10 . C A 0.0 PASS . GT 1/1"
    , "seq 12 . T A 0.0 PASS . GT 1/1"
    , "seq 14 . CGT AGA 0.0 PASS . GT 1/1"
    , "seq 18 . C A 0.0 PASS . GT 1/1"
    , "seq 20 . C A 0.0 PASS . GT 0/1"

  };

  public void testROCTricky() throws IOException, UnindexableDataException {
    checkRoc("tricky", TRICKY_TEMPLATE, EMPTY, FP_TRICKY, FN_TRICKY, false);
  }

  private static final String[] FP_EMPTY = {
      "seq 1 . A T 1 PASS . GT:GQ 1/1:" + PosteriorUtils.phredIfy(1 * MathUtils.LOG_10)
  };

  public void testROCEmpty() throws IOException, UnindexableDataException {
    checkRoc("rocempty", TRICKY_TEMPLATE, EMPTY, FP_EMPTY, EMPTY);
  }
  
  public void testOutsideRef() throws IOException, UnindexableDataException {
    createInput(new String[] {"seq 28 . A G 0.0 PASS . GT 1/1"}, new String[] {"seq 30 . C T 0.0 PASS . GT 0/1"}, new String[] {});
    final File calls = new File(mDir, "calls.vcf.gz");
    final File mutations = new File(mDir, "mutations.vcf.gz");
    final File out = FileUtils.createTempDir("outsideRef", "out", mDir);
    final File template = new File(mDir, "template");

    ReaderTestUtils.getReaderDNA(TEMPLATE, template, null).close();
    final MemoryPrintStream ps = new MemoryPrintStream();
    Diagnostic.setLogStream(ps.printStream());
    try {
      final VcfEvalParams params = VcfEvalParams.builder().baseLineFile(mutations).callsFile(calls)
          .templateFile(template).outputParams(new OutputParams(out, false, false)).create();
      VcfEvalTask.evaluateCalls(params);
    } finally {
      Diagnostic.setLogStream();
    }
    TestUtils.containsAll(ps.toString(),
        "Variant in calls at seq:28 starts outside the length of the reference sequence (27).",
        "Variant in baseline at seq:28 starts outside the length of the reference sequence (27).",
        "Variant in calls at seq:30 starts outside the length of the reference sequence (27).",
        "There were 1 baseline variants skipped due to being too long, overlapping or starting outside the expected reference sequence length.",
        "There were 2 called variants skipped due to being too long, overlapping or starting outside the expected reference sequence length."
        );
  }

  private static final String[] FP_TRICKY_XRX = {
    "seq 1 . A T 1 PASS . GT:GQ 1/1:" + PosteriorUtils.phredIfy(1 * MathUtils.LOG_10)
  , "seq 2 . A AT 9 PASS XRX GT:GQ 1/1:" + PosteriorUtils.phredIfy(9 * MathUtils.LOG_10)
  , "seq 5 . GT G 2 PASS XRX GT:GQ 1/1:" + PosteriorUtils.phredIfy(2 * MathUtils.LOG_10)
  , "seq 10 . CGT AGA 3 PASS XRX GT:GQ 1/1:" + PosteriorUtils.phredIfy(3 * MathUtils.LOG_10)
  , "seq 14 . C A 8 PASS . GT:GQ 1/1:" + PosteriorUtils.phredIfy(8 * MathUtils.LOG_10)
  , "seq 16 . T A 4 PASS . GT:GQ 1/1:" + PosteriorUtils.phredIfy(4 * MathUtils.LOG_10)
  , "seq 20 . C A 10 PASS . GT:GQ 0/1:" + PosteriorUtils.phredIfy(10 * MathUtils.LOG_10)
  };

  private static final String[] FN_TRICKY_XRX = {
    "seq 5 . TC T 0.0 PASS XRX GT 1/1",
    "seq 7 . C CC 0.0 PASS XRX GT 1/1",
    "seq 10 . C A 0.0 PASS . GT 1/1",
    "seq 12 . T A 0.0 PASS . GT 1/1",
    "seq 14 . CGT AGA 0.0 PASS XRX GT 1/1",
    "seq 18 . C A 0.0 PASS . GT 1/1",
    "seq 20 . C A 0.0 PASS . GT 0/1"

  };

  public void testROCTrickyXRX() throws IOException, UnindexableDataException {
    checkRoc("trickyxrx", TRICKY_TEMPLATE, EMPTY, FP_TRICKY_XRX, FN_TRICKY_XRX);
  }

  public void testGetSdfId() {
    Diagnostic.setLogStream();
    assertEquals(new SdfId(0), VcfEvalTask.getSdfId(null));
    final String label = "##TEMPLATE-SDF-ID=";
    final ArrayList<String> header = new ArrayList<>();
    header.add(label + "blahtblah");
    try {
      VcfEvalTask.getSdfId(header);
      fail();
    } catch (final NoTalkbackSlimException e) {
      assertEquals("Invalid header line : " + label + "blahtblah", e.getMessage());
    }
    header.clear();
    header.add(label + "blah");
    try {
      VcfEvalTask.getSdfId(header);
      fail();
    } catch (final NoTalkbackSlimException e) {
      assertEquals("Invalid header line : " + label + "blah", e.getMessage());
    }
    header.clear();
    header.add(label + new SdfId(42).toString());
    assertEquals(new SdfId(42), VcfEvalTask.getSdfId(header));
  }

  public void testPeek() throws IOException {
    Diagnostic.setLogStream();
    BufferedReader reader = new BufferedReader(new StringReader("ooogablah" + StringUtils.LS));
    assertEquals("o", VcfEvalTask.peek(reader));
    reader = new BufferedReader(new StringReader(""));
    assertEquals("", VcfEvalTask.peek(reader));
  }

  public void testCheckHeader() throws IOException {
    final MemoryPrintStream ps = new MemoryPrintStream();
    Diagnostic.setLogStream(ps.printStream());
    try {
      final String label = "##TEMPLATE-SDF-ID=";
      BufferedReader generated = new BufferedReader(new StringReader(""));
      BufferedReader detected = new BufferedReader(new StringReader(""));
      try {
        VcfEvalTask.checkHeader(generated, detected, new SdfId(42));
        fail();
      } catch (final NoTalkbackSlimException ntse) {
        assertEquals("No header found in baseline file", ntse.getMessage());
      }
      ps.reset();
      final SdfId idA = new SdfId();
      SdfId idB;
      do {
        idB = new SdfId();
      } while (idA.check(idB));
      generated = new BufferedReader(new StringReader(label + idA.toString() + StringUtils.LS));
      detected = new BufferedReader(new StringReader(label + idA.toString() + StringUtils.LS));
      VcfEvalTask.checkHeader(generated, detected, idA);
      assertEquals("", ps.toString());

      generated = new BufferedReader(new StringReader(label + idB.toString() + StringUtils.LS));
      detected = new BufferedReader(new StringReader(label + idA.toString() + StringUtils.LS));
      VcfEvalTask.checkHeader(generated, detected, idA);
      assertTrue(ps.toString().contains("Template ID mismatch, baseline variants were not created from the given template"));
      ps.reset();

      generated = new BufferedReader(new StringReader(label + idA.toString() + StringUtils.LS));
      detected = new BufferedReader(new StringReader(label + idB.toString() + StringUtils.LS));
      VcfEvalTask.checkHeader(generated, detected, idA);
      assertTrue(ps.toString().contains("Template ID mismatch, called variants were not created from the given template"));
      ps.reset();

      generated = new BufferedReader(new StringReader(label + idB.toString() + StringUtils.LS));
      detected = new BufferedReader(new StringReader(label + idA.toString() + StringUtils.LS));
      VcfEvalTask.checkHeader(generated, detected, new SdfId(0));
      assertTrue(ps.toString().contains("Template ID mismatch, baseline and called variants were created with different templates"));

      ps.reset();
    } finally {
      Diagnostic.setLogStream();
    }

  }

  public void testGetVariants() throws IOException {
    final File dir = FileUtils.createTempDir("tabixVarianceTest", "test");
    try {
      final File input = new File(dir, "snp_only.vcf.gz");
      FileHelper.resourceToFile("com/rtg/sam/resources/snp_only.vcf.gz", input);
      final File tabix = new File(dir, "snp_only.vcf.gz.tbi");
      FileHelper.resourceToFile("com/rtg/sam/resources/snp_only.vcf.gz.tbi", tabix);
      final File input2 = new File(dir, "snp_only_2.vcf.gz");
      FileHelper.resourceToFile("com/rtg/sam/resources/snp_only.vcf.gz", input2);
      Collection<Pair<String, Integer>> names = new ArrayList<>();
      for (int seq = 1; seq < 32; seq++) {
        names.add(new Pair<>("simulatedSequence" + seq, -1));
      }
      assertTrue(VcfEvalTask.getVariants(VcfEvalParams.builder().callsFile(input).baseLineFile(input).create(), names) instanceof TabixVcfRecordSet);
    } finally {
      FileHelper.deleteAll(dir);
    }
  }

  private static final String TEMPLATE_PREVIOUS_NT = ">chr1" + StringUtils.LS
      + "AAGATGC";

  private static final String[] FN_BUG_PREVIOUS_NT = {
    "chr1 2 . AGA TGA,AATCC 69.9 PASS XRX GT:DP:RE:GQ:AR:AB:RS 0/1:7:0.007:45.5:0.000:0.429:A,3,0.003,C,1,0.001,T,3,0.003"
  };

  private static final String[] TP_BUG_PREVIOUS_NT = {
    "chr1 6 . G C 69.9 PASS XRX GT:DP:RE:GQ:AR:AB:RS 1/1:7:0.007:45.5:0.000:0.429:A,3,0.003,C,1,0.001,T,3,0.003"
  };


  public void testPreviousNt() throws IOException, UnindexableDataException {
    check(TEMPLATE_PREVIOUS_NT, TP_BUG_PREVIOUS_NT, EMPTY, FN_BUG_PREVIOUS_NT, TP_BUG_PREVIOUS_NT, new String[] {}, FN_BUG_PREVIOUS_NT, false, VcfUtils.FORMAT_GENOTYPE_QUALITY);
  }


  private static final String TEMPLATE_COMPLEX_OVERLAP_BUG = ">Chr1" + StringUtils.LS
    + "ACTAAGAAGATTTACAATATGAATTATCCTGATTCTAGTCACCACCACCACAACAACTATCCTCATCTGCTTCTGAGACACTCAAGAAGTTAATAGGTCA" + StringUtils.LS;

  private static final String[] TP_COMPLEX_OVERLAP_BUG = {
    "Chr1 51 . CAACAACTATCCTC CAACAACTATCCTCATCT,CCACAACTATCCTCATCT 1259.6 PASS DP=16;XRX GT:DP:RE:RQ:GQ 1/1:10:0.009:750.7:22.6",
    "Chr1 68 . T TATCT 1192.1 PASS DP=16;XRX GT:DP:RE:RQ:GQ 1/1:9:0.072:753.5:25.6"
  };

  public void testComplexOverlapBug() throws IOException, UnindexableDataException {
    // Test for bug # 1502
    checkRoc("complexoverlap", TEMPLATE_COMPLEX_OVERLAP_BUG, TP_COMPLEX_OVERLAP_BUG, EMPTY, EMPTY);
  }


  private static final String TEMPLATE_BUG = ">chr1" + StringUtils.LS
      + "AAGATG";

  private static final String[] FN_BUG = {
    "chr1 2 . AGA TGA,AATCC 69.9 PASS XRX GT:DP:RE:GQ:AR:AB:RS 0/1:7:0.007:45.5:0.000:0.429:A,3,0.003,C,1,0.001,T,3,0.003"
  };

  private static final String[] TP_BUG = {
    "chr1 6 . G C 69.9 PASS XRX GT:DP:RE:GQ:AR:AB:RS 1/1:7:0.007:45.5:0.000:0.429:A,3,0.003,C,1,0.001,T,3,0.003"
  };

  public void testInfiniteLoop() throws IOException, UnindexableDataException {
    check(TEMPLATE_BUG, TP_BUG, EMPTY, FN_BUG, TP_BUG, new String[] {}, FN_BUG, false, VcfUtils.FORMAT_GENOTYPE_QUALITY);
  }

  /* 123456789012 345
     GTTTTGTTA
      T
      _
      T  TGT
      _  GT_



   */
  private static final String ZERO_COUNT_REF = ">11" + StringUtils.LS
      + "gttttgtta";

  private static final String[] ZERO_COUNT_TP = {
    "11 1 . GT G 53.91 PASS . GT:AD:DP:GQ:PL 0/1:34,7:41:83.90:84,0,910",
  };

  private static final String[] ZERO_COUNT_FP = {
     "11 5 . T G 53.91 PASS . GT:AD:DP:GQ:PL 0/1:34,7:41:83.90:84,0,910",
     "11 6 . G T 58.40 PASS . GT:AD:DP:GQ:PL 0/1:32,9:41:99:132,0,698",
     "11 6 . GT G 255.25 PASS . GT:AD:DP:GQ:PL 0/1:24,12:41:5.19:271,0,5"
  };

  public void testZeroCountBug() throws IOException, UnindexableDataException {
    //the reconstructed reference after applying the FP is same as just applying the TP
    //MutationEval will choose the 3 FP as a TP in comparison to 1 TP as it explains more mutations
    check(ZERO_COUNT_REF, ZERO_COUNT_TP, ZERO_COUNT_FP, EMPTY, ZERO_COUNT_FP, ZERO_COUNT_TP, EMPTY, false, VcfUtils.FORMAT_GENOTYPE_QUALITY);
  }

  static final String REF = ">10"  + StringUtils.LS
      + "CTTTCTCTTTCTTTCTTTCTCTCTTTCCATCCTTCTTTCTTTCTTTCTCTTTCTTTCTTTCTCTCTTTCCATCCTTCTTTCTTTCTTTCCTTCCTTCTTTCTTTCCTTCCTTTCTTTCTTTTTCTTTCTTTCTCTTTCCTTTTTCTTTCTTTCTCTCTCCCTCTTTCTTTTTTTTTTGAGTACCTAAAATCTACTCTCTT"
      + "GGTAAAGTGCCTGCACACAATACAGTATTATTAACTACTCTCCTCATGTTGTACATTGCATCTCTAGACTTCTTCATAGGCTTTGGAAGCCGATAGTTTTGAAAGAACAGCATCATATTATTCAGCATTTCAGAGAGGCCCCCATCAAAAACGAAATGCTCTGTTAGCCTCTGGGCGTGCTGCCAGCTGAGGGCCCCAGC"
      + "TGAGCACCCAGGCCATGGAAGAAAGCTGCCTTGCCCAAGCTCATGTCCCCTCCCGGCGCAGTCCACGCCTGATGTCGGCTTGATGTGGAGAACGTGACCCAGGCTCTGTCTGAATTCAGGACATCCCAGGAGGCATCTCCACGCATCGGAGCTCCCGTGGCCTAGGCTGGAGCCTCCTCGGAGGCCTCTGTGCCACTCCC"
      + "TCCCAGGGTGGGTTATCCTGAGTGCTTCCTAATAGGGGTCTGCAGAGGCAAATCCCAGCCCAGGGAGCCCAACCTCACCCGGGACTTACCCCAGAGAGCATGACAAGGTCACCTGCTCAAAAGGTAGATTCCTAAGCAAAATCAACATTGTGTTCCTGAAGAATGTGGAAAGAAATGTGGGGTAAATCACAGTGCCTGAC"
      + "ATTTTTCTAAAAGGTGATTGGCTTTATCCTTAACATAGATTGTTAAGTCCCACAGCCTTACAAAACAAATTTCAATCACAATCTCAACGGCATTGAGTAAGTTTATATATATATATATATATACAAGTATACATATATATACACACACACACACGTAAAATACACATGGAGTTTTGAAAACTTAGTGCAACAATTTAAAA"
      + "AGGTGAAGTACCTCAATTGCTTTTTCTTATTGACTACTTGTTGAAATGCTATGTTTTTGATAGATTTGGTTAAATAAACATTACCATTATAATTAATCTCCCCTATTTCTCTTTTTTTATTTCTTCATACGGCTAGTAGAACATTTAAAATTGCGTATTTGGCTCTTTTTTACGGCTTGCCTTTTTTTTTCTAGTGAAAG"
      + "ATGGTGTATGATGCTTCTGTCTTTTATGCACTCAACATATTTGTTTCCTAAAATAAGAAAAAAAGAAGTAGCAAAAGGCATAACTCTTCAGAATGAGTCTGATATTAAAAAAAAAAAAAGCCATACATATTTTGTTTCTTCTACTGCTTGAAATTTCAGGTGGCAAATCCCTGGAGAATTTAATTTGTCTTCAATTAGAT"
      + "TTTTAATGTTCCTTTTAGTCACAATCAAAAAGACATTTTCTTAATCAGTGGTGTCTATTACAATAGTTAGGATAATGTAATTTGGCATAATTATATGAAATTAAGATTTAGTACACTCAAGGTTTTTTCCCTCCTTGAATGAAATAAATGCAAAAAGCAAAATAACCCTGAAATCTTACATTTAACAAAAAGTCCCAAAC"
      + "CTATTATTCTTTCCCATTATTGAATTGATGTCTGCAGTGACAGTACTAGGAAGAATAAAATAGCAATCTAGCTTTGAGGCTATAATTTTTAATTGATGGAAGCAGGTACAGGATTTATGCAATATGATATCATCTTAAGAGATGATGGAAGTATAAAACAAGGCAAATGAAACAAACAAATAAACTGTGAAGTGCAGAAA"
      + "GAGAAATTTGAACTGCTTTAACTTCTAACAATCAAGCAGTTGGAATATTTGAGCTGTTACCTGATTATATTTAAGGATTTAGATGTGAATTTAAAAAGTTCCCCAAATAAATCTCCTTGCCATCTCTCAAGTCAGCAGCATTCTGAGAGCTAAGCTGGAAACACACAAAACACATAGCCTGGTTTAATAATAGCTAGACT"
      + "GTGGAAGAAGAAAACTCCTTCACTCTGGCAGGCTAACAAAGTTGATTTAATTTAGCTCTCTTAATTGCTTCATATTGTGAATTTCAACATGCCATGTGTAACAGCTTTCTTTCTAGTGATTTTCTTGTAGAAGTCGTTGTTCCCATTAGGTTTTGAGGTATTTTTCTCTTCATTTCTTATTCTCAAAGGTACTGACTTTA"
      + "TGCTCACTTAAAATTGAATGATTAATGACTCTGTTAGCTACATGATACAATTTTAAAAGAATAAGCCCAAACCACTTAACTAAAATAAATAAGGTTTTTTGTTCATTAACATTAACATAGTATGGTCATAAACTGCCCAATTCCTAAATCACTCTTTATTTCCACTTGGCTTAAAATTATTTTTGCCTATTTAATTTATT"
      + "ATTAGGGAAATTAAACACACACACATTTTTTTCTTGAGCACTCATAATTTGAGATTCATTATGATCACTAAGCAAAAAGTTTACTTATTCAGTCATCATTGCAATGCAACTTAAAAGTCACAACAATCCCAATATAAAAATGTAACTTCATGTAATATCAAAAGACATGTTTTGTTTTATCTTAACAATTATAGTTAATT"
      + "TCCCTGCAGAGTAGATCTCCAGGACTGTCAACAGAATCCCCAAGGACATAGTCATTGGCAGTGGTGGATGCAAAACACACATTAAGATTAATTGCAGAAGCTCTGTCTAAACAGAATGTCTTCTGATATTACTTTATAACTTCTTATATGTAATGCAACTGCCTGAACACTTGCGCTCCTTCATGTGTATAATTCTGTTG"
      + "CACACTGCTTGATTCTGTGTAGGAAACAGTTTGGGAATATTGGGGTTGGGTAAAGTGGAGTCACGCATGTATTACATGCTCATGTTTCCATGTTTCATGCTGATGACGATAAGCCTGCACTGCATGCATCCAATGCATTGCAGATTCACTGGTAAATATATGACCAGTAATCAACACAAATATATTAAATTGCTAGTCAT"
      + "TTTGAATTTTGCTTATGCTGGAATAGGGCTCGTATCCATGAATTATCCAGTTATAACGCAGCTCAGATTCTTGCCGTATGCCTTGTATGTGTGAAACACTTGCCACTTTCAATAATTTTTCACATACAGTTTCCCATTTAACATGTTTGACCAAAAAATAAACACAATATATATATTCTTGACTTAAGACATTGCAACTA"
      + "ATAGAAGCACTCAAATTATGAGAAGATTTCCAATTGTCTATCTACCTAATTCATTCAACAGATAGAAATATTATGCAAGGTAATCTGGTATTTTAGAATGTATAAGAAGCTCTAGTTATTATCTGGATAGAATCAGGAGAATGTACTGATAGGTTTAAATGGATTAAAGACATCTTTCCTTCTTCCTTCAGACCTATAAA"
      + "AAATCAAATTAAATGGAATATAATTGAATGATAAATATAATTTATCTAACATGTCCTTTTGGTCTAGCCTAATATGATAAAGTTCATTGGATTGTTTATGAAATTTCACCTTTTATATTGAGCTTTCCAAGTGTGCTCAGAAAAACAAAAAGGTATATTTTCCAGCCTAAGTTATTAGAAAAATTTAGGACAAATAAAAT"
      + "GAGGTATCATATCTACATCTCTTGCATAACACATTGACCATAGCATGGATTCATGTCTGTTTTCTCTGTTCCAAAATCACAGTGACTATTCAAACCTGAAACTCATTGCTGCTCTCAAATTATAAAAAACCTTTCAAATGGTTCTTACCAATAAATTGTAACCATATTGATCAAGCTAATGAATTTTAACGAAAAAAATT"
      + "GTCTTTCTACATTGTCACATTATAGGTTCCAAAAAGATTCTAGAAATTTTAGATTCCAAAAATTCTAAATGTTTATTTTTCAATTGTTTTCCACTCTATTTGCTCTCCTACTTTTGAGTTTTTAAATTCTTGGGTGGACAGAAGTAAATTAAAATCTCATTTTTATTTATTTATAATTAAATAATCCAAATATATATTTT"
      + "TAGGAAAGCTATTTGATGTTTAGAAAACCTTGTTAGCCTCTTAAACATGTTTTTCTTTGAAGAGACAAGCTGTCTTTTGCCAAGCACAAATTTCAAACACCAACGGCTTCTACGTTCTCTTTGAAACACATAACGTCATTATGATAAATTCTTCTTTGGAGAAGTGAATGCATTTACCCTTCATGTGAAGATGATGCTGC"
      + "AACAGGGTATTTGGCCGTTTTTTTCAGATGGAGAAAGTATCATACATTCCTAATTATTTATCAGAACGATGTGAAATCAAGAAAAAGAAAATAGAAACAGGCGACTGTGGTATTTTTTCAAACATTCAAAGTCAAAAATGTTACCTTATTGTGCGAAATATTCAATTTGACTTCATTTTGACTTCTCTCTCTATTCGTGT"
      + "TCTGTGCCTGCCAGTCTCAAATCCCAAAGCACTGAGTCAGGCTAACGTGAGAAAGTGAAACATATCCACTTTATATTCCAGTTCAATTTTGTTTTAACATGATCTGATGAAAAACTTACGAGACCAAGGTATTAATGAAACTCTAACATTAGTTATTTTTTAATCACCTAATTTTTAGTTATTTTTACTTTTTATTTATT"
      + "CGTTTTTGAGGCATAGTTTATATATTATAAAACTCATTCTTTTAAGCATTCAAATCAGTGGGTTTTAATATATTCACAAAGCGGTGCAATCATCACCACTTTCTAACTCCAGAACATTTTCATAGCCCCAAAAACAAACCACATGCCCTGTGGCACTTCGTCACACCGTTTCTTCTCCTCTCAGGCCCCAGCAATCACTA"
      + "ATCTGTTTCCCGTCTTTATCGATTTGCCCATTCTGGATATTTCATACAGATGAAATAATACAATACACAGCTTTTTGTGTCTGGCTTTCTTCACTTAGCATAACACTTTCAAAGTTTATCCGTATTGTAGCATGTGCCATTACTTCATTTCTTTTTAAAATTGAATAATAATATTTTATTGTGATAAACATATCACATTT"
      + "TGCTTCTCCATGCATCCATTGATGAAAGATTAGTGTCGTGTTCAGTTTTTGGCTTTTTGAATAATGCGCCCATGGATATTCCTGTACAAGCTTCTGAGTGGACCAGTGTTTTCAGTTCTCTTAATATCTACCTGGGAATGGAATTGTTGTAGTCATACAGTAATTTTCCAACTGCCAAGATGTTTCCAAAGTAGCTGCAC"
      + "AATTTTATACTTCCAACAGCAAAGTATGAGGGTTTCAGTTTCTTCATATCCTCATCAACATTTGCAATTGTTTGTCCTTTTTTATTATAAACATGCTAGTGGGTGTCAAGTGGTATCTCATTGTGGTTTTGGTTTGCATTTCCTAATAGTAATCATATTGGCCATGTTTTTATGTACTTAATGGCCATTTGTGTTTATAA"
      + "CTTTCAAAAAGTTTTAGTAATGTTTTTTTCTCATTTTTAATTGGGTTGTTAGCCTTTTTTTTTCTTTTTTTGAGATACAGTCTCACTCTGTCACCTAGGCTGGAGTGCCGTGGTGCGATCTTGGCTCACTGCAACCTCTGCCTTCTGGGTTCAAGCAATTCTTCTGCCTCAGCCTCCCAAGTAGCTGGGACTACAGGTGT"
      + "GTGCCACCAAGCCGAGATAATTTTTGTGGTTTTTTTTAGTAGTGAAGGGGTTTCACCTTGTTGGCCAGGCTGGTTTCGAACTCCTGATCTCAGATGATCCACCCACCTTGGCCTCCCAAAGTGCTGCGATTACAAGCGT"
      ;

  private static final String[] REF_TP = {
    "10 1178 . G T 182.85 PASS . GT:AD:DP:GQ:PL 0/1:11,14:25:99:168,0,223",

  };

  private static final String[] REF_FN = {
    "10 5576 . C T 826.26 PASS . GT:AD:DP:GQ:PL 0/1:37,29:66:99:649,0,933"
  };

  private static final String[] REF_FP = {
    "10 1 . TTTCT T 231.2 PASS . GT:DP:RE:GQ 1/0:29:0.182:55.9",
    "10 31 . C CTTCCTTCCTTCCTTTTTTCTTT 601.9 PASS . GT:DP:RE:GQ 1/0:28:0.188:14.0"

  };

  public void testFun() throws IOException, UnindexableDataException {
    check(REF, REF_TP, REF_FP, REF_FN, REF_TP, REF_FP, REF_FN, false, VcfUtils.FORMAT_GENOTYPE_QUALITY);
  }


  public void testPositionComparator() {
    final VariantPositionComparator vc = new VariantPositionComparator();
    final VcfRecord rec = VcfReader.vcfLineToRecord("chr10 11 . G T 182.85 PASS . GT:AD:DP:GQ:PL 0/1:11,14:25:99:168,0,223".replaceAll(" ", "\t"));
    final DetectedVariant v = new DetectedVariant(rec, 0, RocSortValueExtractor.NULL_EXTRACTOR, false);


    final VcfRecord rec2 = VcfReader.vcfLineToRecord("chr10 12 . G T 182.85 PASS . GT:AD:DP:GQ:PL 0/1:11,14:25:99:168,0,223".replaceAll(" ", "\t"));
    final DetectedVariant v2 = new DetectedVariant(rec2, 0, RocSortValueExtractor.NULL_EXTRACTOR, false);

    final VcfRecord rec3 = VcfReader.vcfLineToRecord("chr10 11 . GT TT 182.85 PASS . GT:AD:DP:GQ:PL 0/1:11,14:25:99:168,0,223".replaceAll(" ", "\t"));
    final DetectedVariant v3 = new DetectedVariant(rec3, 0, RocSortValueExtractor.NULL_EXTRACTOR, false);

    assertEquals(-1, vc.compare(v, v2));
    assertEquals(1, vc.compare(v2, v));

    assertEquals(0, vc.compare(v, v));

    assertEquals(-1, vc.compare(v, v3));
    assertEquals(1, vc.compare(v3, v));
  }
}

