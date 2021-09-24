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
package com.rtg.ngs.longread;


import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.PrintStream;

import com.rtg.index.hash.ngs.OutputProcessor;
import com.rtg.launcher.AbstractNanoTest;
import com.rtg.launcher.MainResult;
import com.rtg.launcher.SequenceParams.SequenceParamsBuilder;
import com.rtg.launcher.globals.CoreGlobalFlags;
import com.rtg.mode.DnaUtils;
import com.rtg.mode.SequenceMode;
import com.rtg.ngs.MapCli;
import com.rtg.ngs.NgsFilterParams;
import com.rtg.ngs.NgsFilterParams.NgsFilterParamsBuilder;
import com.rtg.ngs.NgsOutputParams;
import com.rtg.ngs.NgsParams;
import com.rtg.ngs.NgsParamsBuilder;
import com.rtg.ngs.OutputFilter;
import com.rtg.ngs.PairedEndMapStatistics;
import com.rtg.ngs.SingleEndMapStatistics;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.usage.UsageMetric;
import com.rtg.util.IntegerOrPercentage;
import com.rtg.util.StringUtils;
import com.rtg.util.TestUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.util.io.TestDirectory;

/**
 */
public class LongReadTaskTest extends AbstractNanoTest {

  private static final String READ = "acgtacgtacgtacgtacgtacgtacgtacgtacgacgtacgtacgtacgtacgtacgtacgtacgtacg";
  private static final String TMPL = "acgtacgtacgtacgtacgtacgtacgtacgtacgacgtacgtacgtacgtacgtacgtacgtacgtacgacgtacgtacgtacgtacgtacgtacgtacgtacgacgtacgtacgtacgtacgtacgtacgtacgtacgacgtacgtacgtacgtacgtacgtacgtacgtacgacgtacgtacgtacgtacgtacgtacgtacgtacg";

  protected static NgsParams makeRealParams(File readsl, File readsr, File template, File outputDir, int nt, int wordSize) {
    final SequenceParamsBuilder spb = new SequenceParamsBuilder().directory(readsl).loadNames(false).mode(SequenceMode.UNIDIRECTIONAL).useMemReader(true);

    final NgsParamsBuilder ngsp = NgsParams.builder();

    final IntegerOrPercentage iop = new IntegerOrPercentage(5);

    final SequenceParamsBuilder tspb = new SequenceParamsBuilder().directory(template).mode(SequenceMode.BIDIRECTIONAL);
    final NgsFilterParamsBuilder filter = NgsFilterParams.builder().matedMaxMismatches(iop);
    if (readsr != null) {
      final SequenceParamsBuilder spbr = new SequenceParamsBuilder().directory(readsr).loadNames(false).mode(SequenceMode.UNIDIRECTIONAL).useMemReader(true);
      ngsp.buildSecondParams(spbr.create());
      filter.outputFilter(OutputFilter.PAIRED_END);
    } else {
      filter.outputFilter(OutputFilter.SAM_SINGLE_END);
    }

    final NgsOutputParams outp = NgsOutputParams.builder().keepIntermediateFiles(true).filterParams(filter.create()).outputDir(outputDir).create();
    return ngsp.outputParams(outp).buildFirstParams(spb.create()).searchParams(tspb.create()).useLongReadMapping(true).maxFragmentLength(150).minFragmentLength(5).stepSize(16).numberThreads(nt).create();
  }

  public void testSingleEndMTRev() throws Exception {
    try (final TestDirectory tmpDir = new TestDirectory("longreadtask")) {
      final ByteArrayOutputStream ba = new ByteArrayOutputStream();
      final PrintStream pr = new PrintStream(ba);
      Diagnostic.setLogStream(pr);
      try {
        final File reads = new File(tmpDir, "reads");
        ReaderTestUtils.getReaderDNA(">test\n" + DnaUtils.reverseComplement(READ) + "\n", reads, null);

        final File template = new File(tmpDir, "template");
        ReaderTestUtils.getReaderDNA(">test\n" + TMPL + "\n", template, null);

        final File out = new File(tmpDir, "out");
        assertTrue(out.mkdir());
        try (NgsParams ngsp = makeRealParams(reads, null, template, tmpDir, 2, 16)) {
          final SingleEndMapStatistics stats = new SingleEndMapStatistics(null);
          final UsageMetric usageMetric = new UsageMetric();
          try (OutputProcessor op = ngsp.outputParams().outFilter().makeProcessor(ngsp, stats)) {
            LongReadTask.execute(ngsp.toPositionParams(), op, usageMetric);
            op.finish();
          }

          final String res = TestUtils.stripSAMHeader(FileUtils.fileToString(new File(tmpDir, "alignments.sam")));
          //System.err.println(res);
          mNano.check("longreadtask-testsingleendmtrev.txt", res, false);

          stats.printStatistics(ba);
          pr.flush();

          assertEquals(READ.length(), usageMetric.getMetric());
          TestUtils.containsAll(ba.toString(), "1 100.0% mapped ambiguously (NH > 1)", "Timer total_time", "Index search performance",
            "Timer LR_BS_search", "Timer LR_BS_build_queue", "LR_BS_freeze", "Memory Usage", "LongReadTask Thread create:", "LongReadTask Thread start:", "Timer LR_BS_build_queue_Build-0", "Timer Search_0");
        }
      } finally {
        Diagnostic.setLogStream();
        pr.close();
      }
    }
  }

  public void testSingleEndMTRevWinBits() throws Exception {
    try (final TestDirectory tmpDir = new TestDirectory("longreadtask")) {
      final File reads = new File(tmpDir, "reads");
      ReaderTestUtils.getReaderDNA(">test\n" + DnaUtils.reverseComplement(READ) + "\n", reads, null);

      final File template = new File(tmpDir, "template");
      ReaderTestUtils.getReaderDNA(">test\n" + TMPL + "\n", template, null);

      final File out = new File(tmpDir, "out");
      assertTrue(out.mkdir());
      try (NgsParams ngsp = makeRealParams(reads, null, template, tmpDir, 2, 33)) {
        final SingleEndMapStatistics stats = new SingleEndMapStatistics(null);
        final UsageMetric usageMetric = new UsageMetric();
        try (OutputProcessor op = ngsp.outputParams().outFilter().makeProcessor(ngsp, stats)) {
          LongReadTask.execute(ngsp.toPositionParams(), op, usageMetric);
          op.finish();
        }

        assertEquals(READ.length(), usageMetric.getMetric());
        final String res = TestUtils.stripSAMHeader(FileUtils.fileToString(new File(tmpDir, "alignments.sam")));
        mNano.check("longreadtask-testsingleendmtrev.txt", res, false); // Same result as earlier test
      }
    }
  }

  public void checkSingleEnd(final int threads) throws Exception {
    try (final TestDirectory tmpDir = new TestDirectory("longreadtask")) {
      final ByteArrayOutputStream ba = new ByteArrayOutputStream();
      final PrintStream pr = new PrintStream(ba);
      Diagnostic.setLogStream(pr);
      try {
        final File reads = new File(tmpDir, "reads");
        ReaderTestUtils.getReaderDNA(">test\n" + READ + "\n", reads, null);

        final File template = new File(tmpDir, "template");
        ReaderTestUtils.getReaderDNA(">test\n" + TMPL + "\n", template, null);

        final File out = new File(tmpDir, "out");
        assertTrue(out.mkdir());
        try (NgsParams ngsp = makeRealParams(reads, null, template, tmpDir, threads, 16)) {
          final SingleEndMapStatistics stats = new SingleEndMapStatistics(null);
          final UsageMetric usageMetric = new UsageMetric();
          try (OutputProcessor op = ngsp.outputParams().outFilter().makeProcessor(ngsp, stats)) {
            LongReadTask.execute(ngsp.toPositionParams(), op, usageMetric);
            op.finish();
          }

          final String res = TestUtils.stripSAMHeader(FileUtils.fileToString(new File(tmpDir, "alignments.sam")));
          mNano.check("longreadtask-testsingleend.txt", res, false);

          stats.printStatistics(ba);
          pr.flush();

          //System.err.println(ba.toString());

          assertEquals(READ.length(), usageMetric.getMetric());
          TestUtils.containsAll(ba.toString(), "1 100.0% mapped ambiguously (NH > 1)", "Timer total_time", "Index search performance",
            "Timer LR_BS_search", "Timer LR_BS_build_queue", "LR_BS_freeze", "Memory Usage", "LongReadTask Thread create:", "LongReadTask Thread start:", "Timer LR_BS_build_queue_Build-0", "Timer Search_0");
        }
      } finally {
        Diagnostic.setLogStream();
        pr.close();
      }
    }
  }

  public void testSingleEnd() throws Exception {
    checkSingleEnd(1);
  }

  public void testSingleEndMT() throws Exception {
    checkSingleEnd(2);
  }

  public void testSingleEndMTWinBits() throws Exception {
    try (final TestDirectory tmpDir = new TestDirectory("longreadtask")) {
      final ByteArrayOutputStream ba = new ByteArrayOutputStream();
      final PrintStream pr = new PrintStream(ba);
      Diagnostic.setLogStream(pr);
      try {
        final File reads = new File(tmpDir, "reads");
        ReaderTestUtils.getReaderDNA(">test\n" + READ + "\n", reads, null);

        final File template = new File(tmpDir, "template");
        ReaderTestUtils.getReaderDNA(">test\n" + TMPL + "\n", template, null);

        final File out = new File(tmpDir, "out");
        assertTrue(out.mkdir());
        try (NgsParams ngsp = makeRealParams(reads, null, template, tmpDir, 2, 33)) {
          final SingleEndMapStatistics stats = new SingleEndMapStatistics(null);
          final UsageMetric usageMetric = new UsageMetric();
          try (OutputProcessor op = ngsp.outputParams().outFilter().makeProcessor(ngsp, stats)) {
            LongReadTask.execute(ngsp.toPositionParams(), op, usageMetric);
            op.finish();
          }

          assertEquals(READ.length(), usageMetric.getMetric());
          final String res = TestUtils.stripSAMHeader(FileUtils.fileToString(new File(tmpDir, "alignments.sam")));
          mNano.check("longreadtask-testsingleend.txt", res, false);

        }
      } finally {
        Diagnostic.setLogStream();
        pr.close();
      }
    }
  }

  private static final String LONG_READ_1 = "CATAATGACGGCTGGGCTACTGGACATATGTACGCGGTCTCCGGGAACAAGCAGCATGGAGTTTCCCCTG";
  private static final String LONG_READ_2 = "ACAAGTACAAGGAGGCGGACTATCGCACTCCAAGTTAGCCGGTTTGAACATAGAAGCTCTCCCGAGCTCG";

  private static final String TEMPLATE_LONG = LONG_READ_1 + "atcgag" + LONG_READ_2;

  public void checkPairedEnd(final int threads) throws Exception {
    try (final TestDirectory tmpDir = new TestDirectory("longreadtask")) {
      final File readsl = new File(tmpDir, "readsl");
      ReaderTestUtils.getReaderDNA(">r1\n" + LONG_READ_1 + "\n", readsl, null);
      final File readsr = new File(tmpDir, "readsr");
      ReaderTestUtils.getReaderDNA(">r1\n" + LONG_READ_2 + "\n", readsr, null);

      final File template = new File(tmpDir, "template");
      ReaderTestUtils.getReaderDNA(">test\n" + TEMPLATE_LONG + "\n", template, null);

      final File out = new File(tmpDir, "out");
      assertTrue(out.mkdir());
      try (NgsParams ngsp = makeRealParams(readsl, readsr, template, tmpDir, threads, 16)) {
        final PairedEndMapStatistics stats = new PairedEndMapStatistics(true, null);
        final UsageMetric usageMetric = new UsageMetric();
        try (OutputProcessor op = ngsp.outputParams().outFilter().makeProcessor(ngsp, stats)) {
          LongReadTask.execute(ngsp.toPositionParams(), op, usageMetric);
          op.finish();
        }

        assertEquals(LONG_READ_1.length() + LONG_READ_2.length(), usageMetric.getMetric());
        TestUtils.containsAll(FileUtils.fileToString(new File(tmpDir, "mated.sam")),
          "0\t67\ttest\t1\t55\t70=\t=\t77\t146\tCATAATGACGGCTGGGCTACTGGACATATGTACGCGGTCTCCGGGAACAAGCAGCATGGAGTTTCCCCTG\t*\tAS:i:0\tNM:i:0\tXA:i:0\tIH:i:1",
          "0\t131\ttest\t77\t55\t70=\t=\t1\t-146\tACAAGTACAAGGAGGCGGACTATCGCACTCCAAGTTAGCCGGTTTGAACATAGAAGCTCTCCCGAGCTCG\t*\tAS:i:0\tNM:i:0\tXA:i:0\tIH:i:1");
      }
    }
  }

  public void testPairedEnd() throws Exception {
    checkPairedEnd(1);
  }

  public void testPairedEndMT() throws Exception {
    checkPairedEnd(2);
  }

  private static final String BUG988_TEMPLATE = ">simulatedSequence1" + StringUtils.LS
        + "CGAAGAGAACAGTCCCTCAGCCAGCCTTCTATACAGTGGGCGGGAACAGGCACACGATTTGATGTCTCAATTATGCGCAAATGCATAGATATGAGGCCAAA"
        + "CGACTTCGCCCCCGGATGACTCCGATTGGGAATCAACATGTGGGCCGAGATTGTGTCGTATATGACGATCTGGCCACGGTGTCGGCACTTAATCTACTTGC"
        + "CCTAGTTATATTCCCGACCTAGGTCATCGTTAATTCATCGCTAAGGATCGTGCGGGAATATCATCGGGGTGTCCTGCGAGGCGCCTAATCCCTTAGCCAGA"
        + "AAAATTTGGGTACCTGTTTTAGTAGGTAACTACTGGTTAATTTTCTGTAGTGATTGTTGCTCGTTATCGAAAATGTGCGGCATACTGCTACTAACTCACTT"
        + "ATTTTTGGCGGGTGTTTGTAATAAGATTCTAGGGTTATGTCAGCGATCATACTGCCGTCTTTATTCTCGATATGTGATCAAGGCAACGTCTCATCCAGCAC"
        + "AGATCGACCGTCTCCGTTATATCTCTGTTGACTTTCTTCCTTCTGCTTACTTGATATGCCAACTCCCTAAAGCTTTCGTGGCCTAGAGGCAGAAAGTGTTA"
        + "ATGCGATGACCAGGGTCTGGTAACGATTAGCATCCCTTGTTCGAGTCTCTTAACGATTGCCCTGATGTCGGGGGATCTGAATAGAACCTCCAAACGGTTCT"
        + "AAGTCTGATGGAATTGTGGTTAACTATGTCGGGGTACTCACGTCTGAACCTGCGAAAGGCGGCTTATTAGCGTTCTATTTAGAATGAGTCTGAGGATCTGA"
        + "CTCCCAAAGTTACGTACCTATTCCGAGCACCACGGGCAGTGTCCCCAGGCTGCTTGATTTTAGAGGAGATGCAAACGCCGATAGTGGCAAATTCGCGCTCA"
        + "CCAGAGGTAAAGGGCGGCTGAGCACCCCTTTTGTACTGGGGGTACTCGAAAATCTCGAAATTCCTCTTTCAGGAGGATGGCACGCTTGTAC" + StringUtils.LS;

  private static final String BUG988_READ_1 = "CCAGAGGTAAAGGGCGGCTGAGCACCCCTTTTGTACTGGGGGTACTCGAAAATCTCGAAATTCCTCTTTCAGGAGGATGGCACGCTTGTAC";
  private static final String BUG988_READ_2 =  "CAGAGGTAAAGGGCGGCTGAGCACCCCTTTTGTACTGGGGGTACTCGAAAATCTCGAAATTCCTCTTTCAGGAGGATGGCACGCTTGTACT";
  private static final String BUG988_READ_3 =     "AGGTAAAGGGCGGCTGAGCACCCCTTTTGTACTGGGGGTACTCGAAAATCTCGAAATTCCTCTTTCAGGAGGATGGCACGCTTGTACTATG";
  private static final String BUG988_READ_4 =        "TAAAGGGCGGCTGAGCACCCCTTTTGTACTGGGGGTACTCGAAAATCTCGAAATTCCTCTTTCAGGAGGATGGCACGCTTGTACTCGATGC";
  private static final String BUG988_READ_5 =          "AAGGGCGGCTGAGCACCCCTTTTGTACTGGGGGTACTCGAAAATCTCGAAATTCCTCTTTCAGGAGGATGGCACGCTTGTACTCGATGCGT";
  private static final String BUG988_READ_6 =             "GGCGGCTGAGCACCCCTTTTGTACTGGGGGTACTCGAAAATCTCGAAATTCCTCTTTCAGGAGGATGGCACGCTTGTACTATGGCATTGCA";
  private static final String BUG988_READ_7 =                                 "GTACTGGGGGTACTCGAAAATCTCGAAATTCCTCTTTCAGGAGGATGGCACGCTTGTACTATGCAGTTTAGGCAATTGCAATGGCATGTGG";

  private static final String BUG988_READ_FASTA = ""
          + ">read_1" + StringUtils.LS
          + BUG988_READ_1 + StringUtils.LS
          + ">read_2" + StringUtils.LS
          + BUG988_READ_2 + StringUtils.LS
          + ">read_3" + StringUtils.LS
          + BUG988_READ_3 + StringUtils.LS
          + ">read_4" + StringUtils.LS
          + BUG988_READ_4 + StringUtils.LS
          + ">read_5" + StringUtils.LS
          + BUG988_READ_5 + StringUtils.LS
          + ">read_6" + StringUtils.LS
          + BUG988_READ_6 + StringUtils.LS
          + ">read_7" + StringUtils.LS
          + BUG988_READ_7 + StringUtils.LS;

  public void testBug988() throws IOException {
    try (final TestDirectory dir = new TestDirectory()) {
      final File reference = ReaderTestUtils.getDNADir(BUG988_TEMPLATE, new File(dir, "reference_sdf"));
      final File reads = ReaderTestUtils.getDNADir(BUG988_READ_FASTA, new File(dir, "reads_sdf"));
      final File outdir = new File(dir, "outdir");
      final MemoryPrintStream mps = new MemoryPrintStream();
      final int code = new MapCli().mainInit(new String[] {"-i", reads.getPath(), "-t", reference.getPath(), "-o", outdir.getPath(), "-Z", "--sam", "--read-names", "--unknowns-penalty", "0", "--aligner-mode", "general"}, mps.outputStream(), mps.printStream());
      assertEquals(mps.toString(), 0, code);
      mNano.check("testOffEdge.txt", TestUtils.stripSAMHeader(FileUtils.fileToString(new File(outdir, "alignments.sam"))));
    }
  }

  private static final String OFFEDGE_REF = ""
          + ">simulatedSequence1" + StringUtils.LS
          + "CTTGCTTTGCCTGGCGTGGGACAACCCATATCGGAAATGCAACACTAAGAGGTGTAGATATTACGTCCCATTTACTTCATGGGCTACGTCCGAACCTGTGCAGGTT"
          + "CAACTTTCACCTTACTACCTGTGGCGTTTCCTAATAAACGTGCTCTTACTGGGCGCGCAAAAGTGTCGATGAGCAGATAGTAGTTATTAAAGGCGCTTCCGGAATT"
          + "ATAGTGGGCCCTACCCGGTTGTGGTGAACAAAATTACGACCTGGGAGGTTAACTCTCAATTAAAGATCGTGCGAACACACGGGCACCAAGATGTACGATCGGGAGC"
          + "TCACGTCACATCAGCGATGTCTCATCGTTAAAACAAAGGACGCCAAAAGTTCCTGTGCGTTGATCCTGTTGCGGCGAATACA";

  private static final String OFFEDGE_READS = ""
        + ">read1 60 of beginning 40 off end" + StringUtils.LS
        + "CATCGTTAAAACAAAGGACGCCAAAAGTTCCTGTGCGTTGATCCTGTTGCGGCGAATACACTTGCTTTGCCTGGCGTGGGACAACCCATATCGGAAATGC" + StringUtils.LS
        + ">read2 5 off the end" + StringUtils.LS
        + "TACGATCGGGAGCTCACGTCACATCAGCGATGTCTCATCGTTAAAACAAAGGACGCCAAAAGTTCCTGTGCGTTGATCCTGTTGCGGCGAATACAGGGGG" + StringUtils.LS
        + ">read3 5 off the start" + StringUtils.LS
        + "GGGGGCTTGCTTTGCCTGGCGTGGGACAACCCATATCGGAAATGCAACACTAAGAGGTGTAGATATTACGTCCCATTTACTTCATGGGCTACGTCCGAAC" + StringUtils.LS;

  private static final String OFFEDGE_READS_1 = ""
          + ">read1 5 off front" + StringUtils.LS
          + "GGGGGCTTGCTTTGCCTGGCGTGGGACAACCCATATCGGAAATGCAACACTAAGAGGTGTAGATATTACGTCCCATTTACTTCATGGGCTACGTCCGAAC" + StringUtils.LS;

  private static final String OFFEDGE_READS_2 = ""
          + ">read1 5 off end" + StringUtils.LS
          + "CCCCCTGTATTCGCCGCAACAGGATCAACGCACAGGAACTTTTGGCGTCCTTTGTTTTAACGATGAGACATCGCTGATGTGACGTGAGCTCCCGATCGTA" + StringUtils.LS;

  public void testOffTemplateEndToEnd() throws IOException {
    try (final TestDirectory dir = new TestDirectory()) {
      final File reference = ReaderTestUtils.getDNADir(OFFEDGE_REF, new File(dir, "reference_sdf"));
      final File reads = ReaderTestUtils.getDNADir(OFFEDGE_READS, new File(dir, "reads_sdf"));
      final File outdir = new File(dir, "outdir");
      final MemoryPrintStream mps = new MemoryPrintStream();
      MainResult r = MainResult.run(new MapCli(), "-i", reads.getPath(), "-t", reference.getPath(), "-o", outdir.getPath(), "-Z", "--sam", "--read-names", "--unknowns-penalty", "0", "--aligner-mode", "general", "--soft-clip-distance", "0", "--XX" + CoreGlobalFlags.EDIT_DIST_MISMATCH_SOFT_CLIP, "0");
      assertEquals(r.err(), 0, r.rc());
      mNano.check("testOffEdge2.txt", TestUtils.stripSAMHeader(FileUtils.fileToString(new File(outdir, "alignments.sam"))));
      final File readsPaired = new File(dir, "reads_paired");
      ReaderTestUtils.createPairedReaderDNA(OFFEDGE_READS_1, OFFEDGE_READS_2, readsPaired, null);
      final File outdirPaired = new File(dir, "outdir_paired");
      mps.reset();
      r = MainResult.run(new MapCli(), "-i", readsPaired.getPath(), "-t", reference.getPath(), "-o", outdirPaired.getPath(), "-Z", "--sam", "--read-names", "--unknowns-penalty", "0", "--aligner-mode", "general", "--soft-clip-distance", "0", "--XX" + CoreGlobalFlags.EDIT_DIST_MISMATCH_SOFT_CLIP, "0");
      assertEquals(r.err(), 0, r.rc());
      mNano.check("testOffEdge2paired.txt", TestUtils.stripSAMHeader(FileUtils.fileToString(new File(outdirPaired, "alignments.sam"))));
    }
  }

  /*
  private static final String READ_1 = ">read_l/1" + LS
  + "TGCCAC" + LS
  + ">read_2/2" + LS
  + "TGCTGT" + LS;

  private static final String TEMP_1 = ">template" + LS
  // 1234567890123456789012345678901234
  + "CAGGCAACTGCCACCTTGGTTTTTTGCCCCCCTG" + LS
  //F        111111          1111
  //F       22     22  22    22
  //F        22      22
  //F      11     11           11 11
  //F                           11 11
  //F                            11
  //R11  11            11
  //R      22   22
  //R    22      22
  //R  1111     11
  ;

  private static final String EXP_LONG1TOPN = ""
    + "#read 0 had 12 results with score-indel 0" + LS
    + "template\tR\t1\t1\t0\t0" + LS
    + "template\tF\t1\t6\t0\t0" + LS
    + "template\tF\t1\t13\t0\t0" + LS
    + "template\tF\t1\t30\t0\t0" + LS
    + "template\tF\t1\t33\t0\t0" + LS
    ;

  public void testLong1Topn() {
    final int r = 6, w = 2, a = 0, b = 0, c = 1;
    final boolean cgl = false;
    final int stepSize = 2;
    final boolean useLongReads = true;
    checkLong(new NgsMaskParamsGeneral(r, w, a, b, c, cgl),
        new NgsTestUtils.TestParams(
            useLongReads, stepSize, READ_1, TEMP_1, EXP_LONG1TOPN),
            new NgsTestUtils.NgsFilterPartlyParams(OutputFilter.TOPN, 5), false);
  }
   */
}
