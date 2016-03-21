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
package com.rtg.ngs;

import static com.rtg.util.StringUtils.LS;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.util.Collections;

import com.rtg.launcher.GlobalFlags;
import com.rtg.launcher.SequenceParams;
import com.rtg.mode.SequenceMode;
import com.rtg.reader.PrereadArm;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.reader.SdfId;
import com.rtg.reference.ReferenceGenome;
import com.rtg.usage.UsageMetric;
import com.rtg.util.TestUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.util.io.TestDirectory;
import com.rtg.util.test.FileHelper;
import com.rtg.util.test.NanoRegression;

import junit.framework.TestCase;

/**
 * Tests of Ngs paired end capability.
 */
public class NgsPairedEndTest extends TestCase {

  protected File mDir;

  @Override
  public void setUp() throws Exception {
    Diagnostic.setLogStream();
    mDir = FileHelper.createTempDirectory();
  }

  @Override
  public void tearDown() {
    assertTrue(FileHelper.deleteAll(mDir));
    mDir = null;
  }

  private static final String READ_FIRST = ""
      + ">test" + LS
      + "AGGT" + LS;

  private static final String READ_SECOND = ""
      + ">test" + LS
      + "TAGC" + LS;

  private static final String GENOME = ""
    + ">genome" + LS
    + "AGGT" + "TATT" + "GCTA" + LS
    ;

  //One correctly mated pair - end to end test.
  public void test() throws Exception {
    Diagnostic.setLogStream();
    final int w = 4, a = 0, b = 0, c = 1;
    final boolean cgl = false;
    checkPairedEnd(new NgsMaskParamsGeneral(w, a, b, c),
        new NgsTestUtils.TestPairedEndParams(READ_FIRST, READ_SECOND, GENOME, FileHelper.resourceToString("com/rtg/ngs/resources/ngspe.txt"),
            0, 12, 8L)
        );
  }

  public void testStatistics() throws Exception {
    Diagnostic.setLogStream();
    final int w = 4, a = 0, b = 0, c = 1;
    final boolean cgl = false;
    final NgsMaskParams maskP = new NgsMaskParamsGeneral(w, a, b, c);
    final NgsTestUtils.TestPairedEndParams testP = new NgsTestUtils.TestPairedEndParams(READ_FIRST, READ_SECOND, GENOME, FileHelper.resourceToString("com/rtg/ngs/resources/ngspe.txt"), 0, 12, 8L);
    final ByteArrayOutputStream matedOut = new ByteArrayOutputStream();
    final NgsParams params = getParamsPairedEnd(matedOut, maskP, testP, false);
    final ByteArrayOutputStream stats = new ByteArrayOutputStream();
    final UsageMetric usage = new UsageMetric();
    final NgsTask task = new NgsTask(params, stats, usage);
    try {
      task.exec();
    } finally {
      params.close();
    }
    assertEquals(8, usage.getMetric());
    final MapStatistics statsMap = task.getStatistics();
    assertEquals(2L, statsMap.totalValue(MapStatisticsField.TOTAL_READS));
    assertEquals(2L, statsMap.totalValue(MapStatisticsField.MATED_UNIQUE_READS));
    assertEquals(0L, statsMap.totalValue(MapStatisticsField.MATED_AMBIG_READS));
    assertEquals(0L, statsMap.totalValue(MapStatisticsField.UNMAPPED_NO_HITS));
    assertEquals(0L, statsMap.totalValue(MapStatisticsField.UNMATED_UNIQUE_READS));
    assertEquals(0L, statsMap.value(MapStatisticsField.MISSING, Arm.LEFT));
    assertEquals(0L, statsMap.value(MapStatisticsField.MISSING, Arm.RIGHT));
    assertEquals(100.0, statsMap.totalValueAsPercent(MapStatisticsField.MATED_UNIQUE_READS));
    assertEquals(0.0, statsMap.totalValueAsPercent(MapStatisticsField.MATED_AMBIG_READS));
    assertEquals(0.0, statsMap.totalValueAsPercent(MapStatisticsField.UNMAPPED_NO_HITS));
    assertEquals(0.0, statsMap.totalValueAsPercent(MapStatisticsField.UNMATED_UNIQUE_READS));
  }

  //Test minInsertSize
  public void testMin() throws Exception {
    Diagnostic.setLogStream();
    final int w = 4, a = 0, b = 0, c = 1;
    final boolean cgl = false;
    checkPairedEnd(new NgsMaskParamsGeneral(w, a, b, c),
        new NgsTestUtils.TestPairedEndParams(READ_FIRST, READ_SECOND, GENOME, "", 9, 11, 8L)
        );
  }

  //Test maxInsertSize
  public void testMax() throws Exception {
    Diagnostic.setLogStream();
    final int w = 4, a = 0, b = 0, c = 1;
    final boolean cgl = false;
    checkPairedEnd(new NgsMaskParamsGeneral(w, a, b, c),
        new NgsTestUtils.TestPairedEndParams(READ_FIRST, READ_SECOND, GENOME, "", 0, 8, 8L)
        );
  }

  private static final String GENOME_REVERSE = ""
      + ">genome" + LS
      + "TAGC" + "CATT" + "ACCT" + LS
      ;

  //One correctly mated pair - end to end test - reversed from previous test.
  public void testReverse() throws Exception {
    Diagnostic.setLogStream();
    final int w = 4, a = 0, b = 0, c = 1;
    final boolean cgl = false;
    checkPairedEnd(new NgsMaskParamsGeneral(w, a, b, c),
        new NgsTestUtils.TestPairedEndParams(READ_FIRST, READ_SECOND, GENOME_REVERSE, FileHelper.resourceToString("com/rtg/ngs/resources/ngspe_rev.txt"),
            0,
            12, 8L)
        );
  }

  private static final String GENOME_SIDE = ""
      + ">genome" + LS
      + "ACCT" + "CATT" + "TAGC" + LS
      ;

  //One correctly mated pair - end to end test - sided from previous test.
  public void testSide() throws Exception {
    Diagnostic.setLogStream();
    final int w = 4, a = 0, b = 0, c = 1;
    final boolean cgl = false;
    checkPairedEnd(new NgsMaskParamsGeneral(w, a, b, c),
        new NgsTestUtils.TestPairedEndParams(READ_FIRST, READ_SECOND, GENOME_SIDE, FileHelper.resourceToString("com/rtg/ngs/resources/ngspe_side.txt"),
            0,
            12, 8L)
        );
  }

  private static final String GENOME_FOR_WITHNS = ""
      + ">genome" + LS
      + "AGGT" + "TATT" + "TCAA" + LS
      ;

  private static final String READ_LEFT_WITHNS = ""
      + ">test" + LS
      + "NGGT" + LS;  // the N maps to A.

  private static final String READ_RIGHT_WITHNS = ""
      + ">test" + LS
      + "TNGN" + LS;  // first N maps to A, second N maps to C.

  //One correctly mated pair - end to end test.
  public void testWithNs() throws Exception {
    Diagnostic.setLogStream();
    final int w = 4, a = 0, b = 0, c = 1;
    final boolean cgl = false;
    checkPairedEnd(new NgsMaskParamsGeneral(w, a, b, c),
        new NgsTestUtils.TestPairedEndParams(READ_LEFT_WITHNS, READ_RIGHT_WITHNS, GENOME_FOR_WITHNS, FileHelper.resourceToString("com/rtg/ngs/resources/ngspe_ns.txt"),
            0, 12, 8L)
        );
  }

  // svprep test
  public void testSvPrepFasta() throws Exception {
    final String id = "svprep-expected.txt";
    final String[] additionalArgs = {};
    checkSvPrepFasta(id, additionalArgs);
  }
  // svprep test
  public void testSvPrepSDF() throws Exception {
    final String id = "svprep-expected.txt";
    final String[] additionalArgs = {};
    checkSvPrepSdf(id, additionalArgs);
  }
  // svprep test
  public void testSvPrepNamesFasta() throws Exception {
    final String id = "svprep-expected-names.txt";
    final String[] additionalArgs = {"--read-names"};
    checkSvPrepFasta(id, additionalArgs);
  }
  public void testSvPrepNamesSDF() throws Exception {
    final String id = "svprep-expected-names.txt";
    final String[] additionalArgs = {"--read-names", "-T", "1"};
    checkSvPrepSdf(id, additionalArgs);
  }

  private void checkSvPrepFasta(String id, String[] additionalArgs) throws IOException {
    try (TestDirectory dir = new TestDirectory()) {
      final File leftReads = new File(dir, "left.fasta");
      final File rightReads = new File(dir, "right.fasta");
      final File template = new File(dir, "template");
      final File outDir = new File(dir, "outPutDir");
      ReaderTestUtils.getDNADir(FileHelper.resourceToString("com/rtg/ngs/resources/svprep-template.fasta"), template);
      FileHelper.resourceToFile("com/rtg/ngs/resources/svprep-left.fasta", leftReads);
      FileHelper.resourceToFile("com/rtg/ngs/resources/svprep-right.fasta", rightReads);
      Diagnostic.setLogStream();
      final String[] args = {"-o", outDir.toString()
          , "-l", leftReads.toString()
          , "-r", rightReads.toString()
          , "--format", "fasta"
          , "-t", template.toString()
          , "--sam"
          , "--sam-rg", "@RG\\tPL:ILLUMINA\\tSM:FOO\\tID:BAR"
      };
      checkSvPrepInternal(id, additionalArgs, args, outDir);
    }
  }
  private void checkSvPrepSdf(String id, String[] additionalArgs) throws IOException {
    try (TestDirectory dir = new TestDirectory()) {
      final File reads = new File(dir, "reads");
      ReaderTestUtils.createPairedReaderDNA(FileHelper.resourceToString("com/rtg/ngs/resources/svprep-left.fasta")
          , FileHelper.resourceToString("com/rtg/ngs/resources/svprep-right.fasta")
          , reads, new SdfId());
      final File template = new File(dir, "template");
      final File outDir = new File(dir, "outPutDir");
      ReaderTestUtils.getDNADir(FileHelper.resourceToString("com/rtg/ngs/resources/svprep-template.fasta"), template);
      Diagnostic.setLogStream();
      final String[] args = {"-o", outDir.toString()
          , "-i", reads.toString()
          , "-t", template.toString()
          , "--sam"
          , "--sam-rg", "@RG\\tPL:ILLUMINA\\tSM:FOO\\tID:BAR"
      };
      checkSvPrepInternal(id, additionalArgs, args, outDir);
    }
  }

  private void checkSvPrepInternal(String id, String[] additionalArgs, String[] args, File outDir) throws IOException {
    final NanoRegression r = new NanoRegression(this.getClass());
    try {
      final String[] finalArgs = new String[args.length + additionalArgs.length];
      System.arraycopy(args, 0, finalArgs, 0, args.length);
      System.arraycopy(additionalArgs, 0, finalArgs, args.length, additionalArgs.length);
      final MemoryPrintStream err = new MemoryPrintStream();
      try (ByteArrayOutputStream out = new ByteArrayOutputStream()) {
        GlobalFlags.resetAccessedStatus();
        new MapCli().mainInit(finalArgs, out, err.printStream());
        if (err.toString().length() > 0) {
          System.err.println(err.toString());
        }

        final String full = FileHelper.gzFileToString(new File(outDir, NgsOutputParams.ALIGNMENTS_SAM_FILE_NAME + ".gz"));
        final String actual = TestUtils.stripSAMHeader(full);
        r.check(id, actual);
      }
    } finally {
      r.finish();
    }
  }

  /** the better getParams */
  NgsParams getParamsPairedEnd(final OutputStream out, final NgsMaskParams mask, final NgsTestUtils.TestPairedEndParams test, boolean cg) throws IOException {
    final File subjectsDir = FileHelper.createTempDirectory(mDir);
    final File secondDir = FileHelper.createTempDirectory(mDir);
    final File queriesDir = FileHelper.createTempDirectory(mDir);
    final File hitsDir = FileHelper.createTempDirectory(mDir);
    //System.err.println("hitsDir=" + hitsDir);
    if (cg) {
      ReaderTestUtils.getReaderDNAFastqCG(test.mSubjects, subjectsDir, PrereadArm.LEFT).close();
      ReaderTestUtils.getReaderDNAFastqCG(test.mSubjectsSecond, secondDir, PrereadArm.RIGHT).close();
    } else {
      ReaderTestUtils.getReaderDNA(test.mSubjects, subjectsDir, null).close();
      ReaderTestUtils.getReaderDNA(test.mSubjectsSecond, secondDir, null).close();
    }
    ReaderTestUtils.getReaderDNA(test.mQueries, queriesDir, null).close();
    final SequenceParams subjectParams = SequenceParams.builder().directory(subjectsDir).mode(SequenceMode.UNIDIRECTIONAL).useMemReader(true).create();
    final SequenceParams secondParams = SequenceParams.builder().directory(secondDir).mode(SequenceMode.UNIDIRECTIONAL).useMemReader(true).create();
    final SequenceParams queryParams = SequenceParams.builder().directory(queriesDir).useMemReader(true).loadNames(true).create();
    final NgsFilterParams filterParams = NgsFilterParams.builder().outputFilter(OutputFilter.PAIRED_END).topN(0).errorLimit(10).create();

    final NgsOutputParams outputParams = new NgsTestUtils.OverriddenNgsOutputParams(NgsTestUtils.OverriddenNgsOutputParams.builder().outStream(out).progress(test.mProgress).outputDir(new File(hitsDir, "log")).filterParams(filterParams));
    return NgsParams.builder()
        .useLongReadMapping(false)
        .stepSize(test.mStepSize)
        .buildFirstParams(subjectParams)
        .buildSecondParams(secondParams)
        .searchParams(queryParams)
        .outputParams(outputParams)
        .maskParams(mask)
        .listeners(Collections.singleton(test.mListener))
        .maxFragmentLength(test.mMaxInsertSize)
        .minFragmentLength(test.mMinInsertSize)
        .unknownsPenalty(0)
        .create();
  }

  /** the better check */
  void checkPairedEnd(final NgsMaskParams mask, final NgsTestUtils.TestPairedEndParams test) throws Exception {
    final String actual;
    try (ByteArrayOutputStream out = new ByteArrayOutputStream()) {
      final NgsParams params = getParamsPairedEnd(out, mask, test, false);
      NgsTestUtils.execNgs(params);
      actual = TestUtils.stripSAMHeader(FileUtils.fileToString(new File(params.directory(), NgsOutputParams.MATED_SAM_FILE_NAME)));
    }
    //    System.err.println("Expected\n" + test.mExpected);
    //System.err.println("Actual\n" + actual);
    assertTrue(TestUtils.sameLines(test.mExpected, actual, false));
  }

  private static final String GENOME5 = ""
    + ">genome1" + LS
    + "AGGT" + "TATT" + "GCTA" + LS
    + ">genome2" + LS
    + "AAAA" + "TATT" + "GCTA" + LS
    + ">genome3" + LS
    + "CCCC" + "TATT" + "GCTA" + LS
    + ">genome4" + LS
    + "GGGG" + "TATT" + "GCTA" + LS
    + ">genome5" + LS
    + "TTTT" + "TATT" + "GCTA" + LS
    ;

  private static final String REFERENCE = ""
    + "version\t1" + LS
    + "either\tdef\tdiploid\tlinear" + LS
    + "either\tseq\tgenome1\tdiploid\tlinear" + LS
    + "either\tseq\tgenome2\tdiploid\tlinear" + LS
    + "either\tseq\tgenome3\tdiploid\tlinear" + LS
    + "either\tseq\tgenome4\tdiploid\tlinear" + LS
    + "either\tseq\tgenome5\tdiploid\tlinear" + LS
    ;

  public void testNoSexWithReferenceRelatedBug1596() throws IOException {
    final File testDir = FileUtils.createTempDir("junit", "test");
    try {
      final NanoRegression r = new NanoRegression(this.getClass());
      try {
        final File left = FileHelper.stringToGzFile(READ_FIRST, new File(testDir, "left.fa.gz"));
        final File right = FileHelper.stringToGzFile(READ_SECOND, new File(testDir, "right.fa.gz"));
        try (TestDirectory sdf = new TestDirectory("referencetest1")) {
          ReaderTestUtils.getReaderDNA(GENOME5, sdf, null).close();
          FileUtils.stringToFile(REFERENCE, new File(sdf, ReferenceGenome.REFERENCE_FILE));
          //System.out.println(Arrays.toString(sdf.listFiles()));
          final MemoryPrintStream err = new MemoryPrintStream();
          try (ByteArrayOutputStream out = new ByteArrayOutputStream()) {
            final File outDir = new File(testDir, "out");
            GlobalFlags.resetAccessedStatus();
            new MapCli().mainInit(new String[] {"-t", sdf.getPath(), "-l", left.getPath(), "-r", right.getPath(),
              "--format", "fasta", "-w", "4", "-c", "1", "-a", "0", "-b", "0", "--sam",
              "--sam-rg", "@RG\\tID:READGROUP1\\tSM:BACT_SAMPLE\\tPL:ILLUMINA", "-o", outDir.getPath()}, out, err.printStream());
            if (err.toString().length() > 0) {
              System.err.println(err.toString());
            }
            //System.out.println(Arrays.toString(outDir.listFiles()));
            final String log = FileUtils.fileToString(new File(outDir, "map.log"));
            assertFalse(log, log.contains("tarting check"));
            final String full = FileHelper.gzFileToString(new File(outDir, NgsOutputParams.ALIGNMENTS_SAM_FILE_NAME + ".gz"));
            final String actual = TestUtils.stripSAMHeader(full);
            r.check("no-sex-ref.sam", actual);
          }
        }
      } finally {
        r.finish();
      }
    } finally {
      FileHelper.deleteAll(testDir);
    }
  }
}
