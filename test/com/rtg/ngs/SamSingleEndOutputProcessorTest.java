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


import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.PrintStream;

import com.rtg.AbstractTest;
import com.rtg.bed.BedUtils;
import com.rtg.index.hash.ngs.OutputProcessor;
import com.rtg.launcher.HashingRegion;
import com.rtg.launcher.MockSequenceParams;
import com.rtg.launcher.SequenceParams;
import com.rtg.mode.SequenceMode;
import com.rtg.mode.SequenceType;
import com.rtg.reader.MockSequencesReader;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.util.IORunnable;
import com.rtg.util.SimpleThreadPool;
import com.rtg.util.StringUtils;
import com.rtg.util.TestUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.IOUtils;
import com.rtg.util.io.TestDirectory;
import com.rtg.util.test.FileHelper;

/**
 * Tests corresponding class
 */
public class SamSingleEndOutputProcessorTest extends AbstractTest {

  protected File mDir = null;
  @Override
  public void setUp() throws IOException {
    mDir = FileHelper.createTempDirectory();
    super.setUp();
  }

  @Override
  public void tearDown() throws IOException {
    assertTrue(mDir == null || !mDir.exists() || FileHelper.deleteAll(mDir));
    mDir = null;
    super.tearDown();
  }

  static final String TEMPLATE = ">t" + StringUtils.LS + "tgcaagacaagagggcctcc" + StringUtils.LS;
  static final String TEMP1 = "tgcaagacaagagggcctcc";
  static final String TEMP2 = "ggaggccctcttgtcttgca";
  static final String READS = ">r1" + StringUtils.LS + TEMP1 + StringUtils.LS
  + ">r2" + StringUtils.LS + TEMP2 + StringUtils.LS;

  protected NgsParamsBuilder getDefaultBuilder(String reads) throws IOException {
    return getDefaultBuilder(null, false, reads);
  }

  private NgsParamsBuilder getDefaultBuilder(final File tempDir, final boolean gzipOutputs, String reads) throws IOException {
    final File templateok = FileUtils.createTempDir("template", "ngs", mDir);
    final File readsok = FileUtils.createTempDir("reads", "ngs", mDir);
    final File hitsDir = new File(mDir, "hitDir");

    ReaderTestUtils.getReaderDNA(TEMPLATE, templateok, null).close();
    ReaderTestUtils.getReaderDNA(reads, readsok, null).close();

    final NgsFilterParams filterParams = NgsFilterParams.builder().outputFilter(OutputFilter.PAIRED_END)
    .topN(10).errorLimit(5).zip(gzipOutputs).create();
    final NgsOutputParams outputParams = NgsOutputParams.builder()
    .progress(false).outputDir(hitsDir)
    .tempFilesDir(tempDir).filterParams(filterParams).create();

    return NgsParams.builder()
    .buildFirstParams(SequenceParams.builder().directory(readsok).useMemReader(true).create())
    .buildSecondParams(SequenceParams.builder().directory(readsok).useMemReader(true).create())
    .searchParams(SequenceParams.builder().directory(templateok).useMemReader(true).loadNames(true).create())
    .outputParams(outputParams)
    .substitutionPenalty(1).gapOpenPenalty(1).gapExtendPenalty(1).unknownsPenalty(1)
    .maxFragmentLength(1000).minFragmentLength(0);
  }

  //  static class DummyPairedEndOutputProcessor extends PairedEndOutputProcessor {
  //
  //    public DummyPairedEndOutputProcessor(final NgsParams param) {
  //      super(param, null, param.outStream());
  //    }
  //  }

  public void test() throws Exception {
    try (NgsParams param = getDefaultBuilder(READS).create()) {
      final ByteArrayOutputStream log = new ByteArrayOutputStream();
      final PrintStream prLog = new PrintStream(log);
      Diagnostic.setLogStream(prLog);
      try {
        final SamSingleEndOutputProcessor sseop = new SamSingleEndOutputProcessor(param, null, false);
        TestUtils.containsAll(sseop.toString(), "SamSingleEndOutputProcessor; topn= TopNImplementation");
        assertTrue(sseop.toString().contains("TopNImplementationSync"));
        final OutputProcessor cloned = sseop.threadClone(HashingRegion.NONE);
        cloned.process(0, "F", 0, 1, 0, 0);
        cloned.threadFinish();
        //sseop.process(0, "R", 1, 16, 0, 0);

        sseop.finish();
        sseop.close();
        ///assertNull(peop.mSamWriter);
        final String outStr = IOUtils.readAll(new File(param.outputParams().directory(), "alignments.sam"));

        assertTrue(outStr.contains("0\t0\tt\t1\t37\t20=\t*\t0\t0\tTGCAAGACAAGAGGGCCTCC\t*\tAS:i:0\tNM:i:0\tIH:i:1\tNH:i:1"));
        //assertTrue(outStr.contains("0\t99\tt\t1\t255\t20=\t=\t1\t20\tTGCAAGACAAGAGGGCCTCC\t*\tAS:i:0\tNM:i:0\tMQ:i:255"));

        sseop.mUnmappedTracker.calculateStatistics(false, false);
        prLog.flush();

        //        TestUtils.containsAll(baos.toString(), new String[] {
        //          "Sliding window collector statistics",
        //          "Setting stats:",
        //          "hits = 2"
        //        });

      } finally {
        prLog.close();
        //System.err.println("log:\n" + log.toString());

        Diagnostic.setLogStream();
      }
      TestUtils.containsAll(log.toString(), "Extracting hits for reads",
        "AlignmentOutput",
        "Processing ",
        " hits for reads",
        "Merging alignment results");
    }
  }

  static final int MAX_COORD = 50;

  private static class SimpleProcess implements IORunnable {
    private final OutputProcessor mParentProc;
    private final HashingRegion mRegion;
    private OutputProcessor mProc;
    private final int mThreadNum;

    SimpleProcess(final OutputProcessor proc, HashingRegion region, final int threadNum) {
      mParentProc = proc;
      mThreadNum = threadNum;
      mRegion = region;
    }

    @Override
    public void run() {
      try {
        mProc = mParentProc.threadClone(mRegion);
        try {
          mProc.process(0, "F", 0, 1, mThreadNum, 0);
          mProc.process(0, "R", 1, 5, mThreadNum, 0);
          mProc.process(0, "F", 0, 2, mThreadNum, 0);
          mProc.process(0, "R", 1, 8, mThreadNum, 0);
          mProc.process(0, "F", 0, 2, mThreadNum, 0);
          mProc.process(0, "R", 1, 29, mThreadNum, 0);
          mProc.process(0, "F", 0, 19, mThreadNum, 0);
          mProc.process(0, "R", 1, 39, mThreadNum, 0);
        } finally {
          mProc.threadFinish();
        }
      } catch (final IOException e) {
        fail(e.getMessage());
      }
    }
  }

  public void testThreadClone1() throws Exception {
    checkThreadClone(null, false);
  }

  public void testThreadClone2() throws Exception {
    checkThreadClone(null, true);
  }
  //
  //  public void testThreadClone3() {
  //    final File tmp = FileUtils.createTempDir("topnsync", "clone");
  //    try {
  //      checkThreadClone(tmp, true);
  //    } finally {
  //      assertTrue(FileUtils.deleteAll(tmp));
  //    }
  //  }

  //  public void testThreadCloneLocalZipText() {
  //    final File tmp = FileUtils.createTempDir("topnsync", "clone");
  //    try {
  //      checkThreadClone(tmp, false);
  //    } finally {
  //      assertTrue(FileUtils.deleteAll(tmp));
  //    }
  //  }

  static final String TB = "\t";
  static final String EXPECTED = ""
    + "0" + TB + "0" + TB + "t" + TB + "1" + TB + "37" + TB + "20=" + TB + "*" + TB + "0" + TB + "0" + TB + "TGCAAGACAAGAGGGCCTCC" + TB + "*" + TB + "AS:i:0" + TB + "NM:i:0" + TB + "IH:i:1\tNH:i:1" + StringUtils.LS
    + "1" + TB + "16" + TB + "t" + TB + "1" + TB + "37" + TB + "20=" + TB + "*" + TB + "0" + TB + "0" + TB + "TGCAAGACAAGAGGGCCTCC" + TB + "*" + TB + "AS:i:0" + TB + "NM:i:0" + TB + "IH:i:1\tNH:i:1" + StringUtils.LS;

  /**
   * Copied from TopnPairedEndOutputProcessorSyncTest of threadClone method,
   *
   * @param tempDir null or the name of a directory for temporary SAM files
   * @param gzipResults true means compress the final result SAM files
   * @throws Exception on IO error
   */
  public void checkThreadClone(final File tempDir, final boolean gzipResults) throws Exception {
    final int numThreads = 4;
    final SimpleThreadPool stp = new SimpleThreadPool(numThreads, "TestSamSingleEnd", true);
    final NgsParams params = getDefaultBuilder(tempDir, gzipResults, READS).numberThreads(numThreads).create();
    final ByteArrayOutputStream log = new ByteArrayOutputStream();
    final PrintStream prLog = new PrintStream(log);
    Diagnostic.setLogStream(prLog);

    final SamSingleEndOutputProcessor sync = new SamSingleEndOutputProcessor(params, null, true);
    TestUtils.containsAll(sync.toString(), "SamSingleEndOutputProcessor; topn= TopNImplementationSync",
    "temp files gzipped= ");
    //todo True / true
    try {
      for (int i = 0; i < numThreads; i++) {
        final long padding = params.calculateThreadPadding();
        final long start = i * MAX_COORD / numThreads;
        final long end = (i + 1) * MAX_COORD / numThreads;
        final HashingRegion region = new HashingRegion(0, start, 0, end, Math.max(0, start - padding), Math.min(MAX_COORD, end + padding));
        //System.err.println(region);
        stp.execute(new SimpleProcess(sync, region, i));
      }
      stp.terminate();
      //System.err.println(sync.toString());
      TestUtils.containsAll(sync.toString(), "SamSingleEndOutputProcessor; topn= TopNImplementationSync",
          "numsequences= 2", "read= 0 stats= C", "read= 1 stats= C");
      sync.finish();
    } finally {
      sync.close();
      params.close();
      prLog.flush();
      prLog.close();
      //System.err.println(log.toString());
      if (tempDir == null) {
        TestUtils.containsAll(log.toString(), "TEMP_SAM_0.gz", "TEMP_SAM_1.gz", "TEMP_SAM_2.gz", "TEMP_SAM_3.gz",
        "Writing unmapped records");
      }
      Diagnostic.setLogStream();
    }
    final File gzFile = new File(new File(mDir, "hitDir"), NgsOutputParams.ALIGNMENTS_SAM_FILE_NAME + FileUtils.GZ_SUFFIX);
    final File outFile = new File(new File(mDir, "hitDir"), NgsOutputParams.ALIGNMENTS_SAM_FILE_NAME);
    if (gzipResults) {
      assertTrue(gzFile.exists());
      assertTrue(!outFile.exists());
    } else {
      assertTrue(outFile.exists());
      assertTrue(!gzFile.exists());

    }
    final String contents = gzipResults
    ? TestUtils.stripSAMHeader(FileHelper.gzFileToString(new File(new File(mDir, "hitDir"), NgsOutputParams.ALIGNMENTS_SAM_FILE_NAME + FileUtils.GZ_SUFFIX)))
        : TestUtils.stripSAMHeader(FileUtils.fileToString(new File(new File(mDir, "hitDir"), NgsOutputParams.ALIGNMENTS_SAM_FILE_NAME)));
    assertTrue(TestUtils.sameLines(EXPECTED, contents, false));
  }

  public void testCalibrateUnmappedFile() throws IOException {
    try (final TestDirectory tempDir = new TestDirectory()) {
      final NgsParams dummy = NgsParams.builder().outputParams(NgsOutputParams.builder().calibrate(true).create()).create();
      final File unmappedFile = new File(tempDir, "unmapped.sam.gz");
      final File unmappedCalFile = new File(tempDir, "unmapped.sam.gz.calibration");
      SamSingleEndOutputProcessor.calibrateUnmappedFile(dummy, unmappedFile);
      assertTrue(unmappedCalFile.isFile());
      final String output = FileUtils.fileToString(unmappedCalFile);
      assertEquals(3, output.split("\n").length);
      TestUtils.containsAll(output, "#Version", "#CL", "@covar\treadgroup\tbasequality\tsequence\tequal\tdiff\tins\tdel");
    }
  }

  public void testCalibrateUnmappedFileWithBedRegions() throws IOException {
    try (final TestDirectory tempDir = new TestDirectory()) {
      final File regionsBedFile = FileUtils.stringToFile("seq0\t0\t5", new File(tempDir, "regions.bed"));
      final NgsParams dummy = NgsParams.builder()
      .searchParams(new MockSequenceParams(new MockSequencesReader(SequenceType.DNA, 2, 10), SequenceMode.BIDIRECTIONAL))
      .outputParams(NgsOutputParams.builder().calibrate(true).calibrateRegions(BedUtils.regions(regionsBedFile)).create()).create();

      final File unmappedFile = new File(tempDir, "unmapped.sam.gz");
      final File unmappedCalFile = new File(tempDir, "unmapped.sam.gz.calibration");
      SamSingleEndOutputProcessor.calibrateUnmappedFile(dummy, unmappedFile);
      assertTrue(unmappedCalFile.isFile());
      final String output = FileUtils.fileToString(unmappedCalFile);
      assertEquals(5, output.split("\n").length);
      TestUtils.containsAll(output, "#Version", "#CL", "@sequence\t5\tseq0", "@sequence\t0\tseq1", "@covar\treadgroup\tbasequality\tsequence\tequal\tdiff\tins\tdel");
    }
  }

}
