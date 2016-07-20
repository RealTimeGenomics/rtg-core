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
import java.io.OutputStream;
import java.io.PrintStream;

import com.rtg.index.hash.ngs.OutputProcessor;
import com.rtg.launcher.HashingRegion;
import com.rtg.launcher.SequenceParams;
import com.rtg.launcher.globals.GlobalFlags;
import com.rtg.ngs.blocking.MapQScoringReadBlocker;
import com.rtg.ngs.tempstage.AbstractTempFileWriter;
import com.rtg.ngs.tempstage.SingleEndTempFileWriter;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.reader.SequencesReader;
import com.rtg.util.IORunnable;
import com.rtg.util.SimpleThreadPool;
import com.rtg.util.StringUtils;
import com.rtg.util.TestUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.FileUtils;
import com.rtg.util.test.FileHelper;
import com.rtg.util.test.NanoRegression;

import junit.framework.TestCase;

/**
 */
public class DummyAlignmentWriterThreadTest extends TestCase { // PairedEndOutputProcessorTest {

  protected File mDir = null;
  protected NanoRegression mNano = null;

  @Override
  public void setUp() throws Exception {
    mDir = FileHelper.createTempDirectory();
    mNano = new NanoRegression(this.getClass());
    Diagnostic.setLogStream();
  }

  @Override
  public void tearDown() throws Exception {
    GlobalFlags.resetAccessedStatus();
    assertTrue(mDir == null || !mDir.exists() || FileHelper.deleteAll(mDir));
    try {
      mNano.finish();
    } finally {
      mNano = null;
    }
    mDir = null;
  }

  private static final class MyAlignmentWriterThread extends AbstractAlignmentWriterThread {
    MyAlignmentWriterThread(AbstractTempFileWriter writer, MatchResult results, long chunkStart, long chunkEnd, HashingRegion region, int threadNumber) {
      super(writer, results, chunkStart, chunkEnd, region, threadNumber);
      //System.err.println(this.makeString());
      assertEquals(this.getName(), "Alignment Processing Thread ");

    }

    String mResultString = "";
    @Override
    protected void handleResult(int templateId, int encodedReadId, int position, boolean reverse) {
      //final long readId = result.mEncodedReadId;
      mResultString += "templateId=" + templateId
      + " position=" + position
      + " readId=" + encodedReadId
      + " reverse=" + reverse + LS;

      //((SamSingleEndAlignmentWriter) mSamWriter).alignmentResult(readId, result.mReverse, result.mPosition);
    }

    public String stringResults() {
      return mResultString;
    }
  }

  static final OutputStream NULL_OUTPUTSTREAM = TestUtils.getNullOutputStream();

  public void testMyAlignmentWriterThread() throws IOException {
    final File tempDir = FileUtils.createTempDir("dummyalignmentwriterthread", "unmatedlocal");
    try {
      final int numThreads = 4;
      NgsParams params = getDefaultBuilder(tempDir, false, TEMPLATE2).numberThreads(numThreads).create();
      final SamSingleEndOutputProcessor op = new SamSingleEndOutputProcessor(params, null, true);
      final SingleEndTempFileWriter writer = new SingleEndTempFileWriter(params, op.mUnmappedTracker, SharedResources.generateSharedResources(params));
      final MapQScoringReadBlocker blockerLeft = new MapQScoringReadBlocker((int) params.buildFirstParams().reader().numberSequences(), 10);
      writer.initialiseAlignments(NULL_OUTPUTSTREAM, blockerLeft);
      final MatchResult results = new MatchResult(0);
      final MyAlignmentWriterThread awt = new MyAlignmentWriterThread(writer, results, -1, -1, HashingRegion.NONE, 1);
      assertEquals("1", awt.toString());

      TestUtils.containsAll(awt.makeString(), "results.length= 0",
        "thread number= 1",
        "chunk start= -1",
        "chunk end= -1",
        "clip end position= -1",
        "end id= -1",
        "start id= -1",
        "clip start position= -1");

      params = getDefaultBuilder(tempDir, false, TEMPLATE2).numberThreads(numThreads).create();
      MatchResult results3 = new MatchResult(1);
      results3.addMatchResult(0, 5, 0, false);
      //assertEquals("templateId=0 position=0 readId=0 reverse=" + Boolean.FALSE, results3[0].toString());

      final MyAlignmentWriterThread awt3 = new MyAlignmentWriterThread(writer, results3, 1, 1, HashingRegion.NONE, 3);
      TestUtils.containsAll(awt3.makeString(), "clip start position= -1",
        "chunk start= 1",
        "chunk end= 1",
        "clip end position= -1",
        "end id= -1",
        "start id= -1");

      results3 = new MatchResult(0);
      results3.addMatchResult(33, 5, 0, false);
      //assertEquals("templateId=0 position=0 readId=0 reverse=" + Boolean.FALSE, results3[0].toString());
      final MyAlignmentWriterThread awt4 = new MyAlignmentWriterThread(writer, results3, 0, 1, new HashingRegion(33, 0, 33, 1065, 0, 0), 0);

      //System.err.println(awt4.makeString());
      TestUtils.containsAll(awt4.makeString(), "clip start position= 0",
        "chunk start= 0",
        "chunk end= 1",
        "clip end position= 1065",
        "end id= 33",
        "start id= 33");

      final MatchResult results4 = new MatchResult(20);
      for (int i = 0; i < 10; i++) {
        results4.addMatchResult(0, i, 0, false);
      }
      for (int i = 10; i < 20; i++) {
        results4.addMatchResult(1, i, 0, false);
      }
      //assertEquals("templateId=0 position=0 readId=0 reverse=" + Boolean.FALSE, results3[0].toString());
      final MyAlignmentWriterThread awt5 = new MyAlignmentWriterThread(writer, results4, 0, 20, HashingRegion.NONE, 1);
      TestUtils.containsAll(awt5.makeString(), "chunk end= 20",
        "chunk start= 0");
      awt5.run();
      assertEquals(20, awt5.stringResults().split(StringUtils.LS).length);

      final SingleEndTempFileWriter writer6 = new SingleEndTempFileWriter(params, op.mUnmappedTracker, SharedResources.generateSharedResources(params));
      final MyAlignmentWriterThread awt6 = new MyAlignmentWriterThread(writer6, results4, 10, 20, HashingRegion.NONE, 2);
      TestUtils.containsAll(awt6.makeString(), "chunk end= 20",
        "chunk start= 10");
      final ByteArrayOutputStream log = new ByteArrayOutputStream();
      final PrintStream prLog = new PrintStream(log);

      Diagnostic.setLogStream(prLog);
      awt6.run();
      prLog.close();
      TestUtils.containsAll(log.toString(), "starting", "finished");
      Diagnostic.setLogStream();
      assertEquals(10, awt6.stringResults().split(StringUtils.LS).length);
    } finally {
      assertTrue(FileHelper.deleteAll(tempDir));
    }
  }

  private static class SimpleProcess2 implements IORunnable {

    private final OutputProcessor mParentProc;
    private final HashingRegion mRegion;
    private OutputProcessor mProc;
    private final int mThreadNum;
    SimpleProcess2(final OutputProcessor proc, HashingRegion region, final int threadNum) {
      mParentProc = proc;
      mRegion = region;
      mThreadNum = threadNum;
    }

    @Override
    public void run() {
      try {
        mProc = mParentProc.threadClone(mRegion);
        try {
          mProc.process(0, "F", 0, 1, mThreadNum, 0);
        } finally {
          mProc.threadFinish();
        }
      } catch (final IOException e) {
        fail();
      }
    }
  }

  public void testUnmatedLocal() throws Exception {
    final File tmp = FileUtils.createTempDir("dummyalignmentwriterthread", "unmatedlocal");
    try {
      checkUnmated(tmp);
    } finally {
      assertTrue(FileHelper.deleteAll(tmp));
    }
  }

  public void testUnmatedRemote() throws Exception {
    checkUnmated(null);
  }

  /**
   * Test of threadClone method, of class PairedEndOutputProcessorSync.
   */
  public void checkUnmated(final File tempDir) throws Exception {
    final ByteArrayOutputStream logBytes = new ByteArrayOutputStream();
    try (PrintStream log = new PrintStream(logBytes)) {
      Diagnostic.setLogStream(log);
      try {
        final int numThreads = 4;
        final SimpleThreadPool stp = new SimpleThreadPool(numThreads, "DummyUnmated", true);
        final NgsParams params = getDefaultBuilder(tempDir, false, TEMPLATE).numberThreads(numThreads).create();
        final SequencesReader ref = params.searchParams().reader();
        final HashingRegion[] regions = HashingRegion.splitWorkload(ref, params.sex(), 0, ref.numberSequences(), params.numberThreads() * params.threadMultiplier(), HashingRegion.DEFAULT_MIN_CHUNK_SIZE, params.calculateThreadPadding());
        try (TopNPairedEndOutputProcessorSync sync = new TopNPairedEndOutputProcessorSync(params, null, true, true)) {
          for (int i = 0; i < regions.length; i++) {
            stp.execute(new SimpleProcess2(sync, regions[i], i));
          }
          stp.terminate();
          sync.finish();
        }
        final String contentsUnmapped = TestUtils.stripSAMHeader(FileUtils.fileToString(new File(new File(mDir, "hitDir"), NgsOutputParams.UNMAPPED_SAM_FILE_NAME)));
        final String contentsUnmated = TestUtils.stripSAMHeader(FileUtils.fileToString(new File(new File(mDir, "hitDir"), NgsOutputParams.UNMATED_SAM_FILE_NAME)));
        mNano.check("dawtt-unmated", contentsUnmated, false);
        mNano.check("dawtt-unmapped", contentsUnmapped, false);
      } finally {
        Diagnostic.setLogStream();
      }
    }
  }

  static final String TEMPLATE = ">t" + StringUtils.LS + "tgcaagacaagagggcctcc" + StringUtils.LS;
  static final String TEMPLATE2 = ">t0" + StringUtils.LS + "tgcaagacaagagggcctcc" + StringUtils.LS
  + ">t1" + StringUtils.LS + "tgcaagacaagagggcctcc" + StringUtils.LS;
  static final String TEMP_LEFT = "tgcaagacaagagggcctcc";
  static final String TEMP_RIGHT = "ggaggccctcttgtcttgca";
  static final String READ_LEFT = ">r" + StringUtils.LS + TEMP_LEFT + StringUtils.LS;
  static final String READ_RIGHT = ">r" + StringUtils.LS + TEMP_RIGHT + StringUtils.LS;

  private NgsParamsBuilder getDefaultBuilder(final File tempDir, final boolean gzipOutputs, String template) throws IOException {
    final File templateok = FileUtils.createTempDir("template", "ngs", mDir);
    final File leftok = FileUtils.createTempDir("left", "ngs", mDir);
    final File rightok = FileUtils.createTempDir("right", "ngs", mDir);
    final File hitsDir = new File(mDir, "hitDir");

    ReaderTestUtils.getReaderDNA(template, templateok, null).close();
    ReaderTestUtils.getReaderDNA(READ_LEFT, leftok, null).close();
    ReaderTestUtils.getReaderDNA(READ_RIGHT, rightok, null).close();

    final NgsFilterParams filterParams = NgsFilterParams.builder().outputFilter(OutputFilter.TOPN_PAIRED_END)
    .topN(10).errorLimit(5).zip(gzipOutputs).create();

    //
    final NgsMaskParams maskParams = new NgsMaskParamsGeneral(4, 1, 2, 2);
    final NgsOutputParams outputParams = NgsOutputParams.builder()
    .progress(false).outputDir(hitsDir)
    .tempFilesDir(tempDir).filterParams(filterParams).create();

    return NgsParams.builder()
    .buildFirstParams(SequenceParams.builder().directory(leftok).useMemReader(true).create())
    .buildSecondParams(SequenceParams.builder().directory(rightok).useMemReader(true).create())
    .searchParams(SequenceParams.builder().directory(templateok).useMemReader(true).loadNames(true).create())
    .outputParams(outputParams)
    .maskParams(maskParams)
    .substitutionPenalty(1).gapOpenPenalty(1).gapExtendPenalty(1).unknownsPenalty(0)
    .maxFragmentLength(1000).minFragmentLength(0);
  }

  private static final MatchResult RESULTS = initResults1();
  static MatchResult initResults1() {
    final MatchResult ret = new MatchResult(10);
    ret.addMatchResult(0, 1, 0, false);
    ret.addMatchResult(1, 1, 1, false);
    ret.addMatchResult(2, 1, 2, false);
    ret.addMatchResult(3, 1, 3, false);
    ret.addMatchResult(4, 1, 4, false);
    ret.addMatchResult(5, 1, 5, false);
    ret.addMatchResult(6, 1, 6, false);
    ret.addMatchResult(7, 1, 7, false);
    ret.addMatchResult(8, 1, 8, false);
    ret.addMatchResult(9, 1, 9, false);
    return ret;
  }

  public void testRegions() {
    assertEquals(new HashingRegion(7, 1, 9, 34, -1, -1), getRegion(8, 33, RESULTS, 7));
    assertEquals(new HashingRegion(0, 0, 1, 1, -1, -1), getRegion(8, 33, RESULTS, 0));
    assertEquals(new HashingRegion(5, 1, 6, 1, -1, -1), getRegion(8, 33, RESULTS, 5));
  }

  public void testThreadsGTResults() {
    assertEquals(new HashingRegion(0, 0, 1, 1, 0, 34), getRegion(100, 33, RESULTS, 0));
    for (int i = 1; i < 9; i++) {
      assertEquals(new HashingRegion(i, 1, i + 1, 1, 0, 34), getRegion(100, 33, RESULTS, i));
    }
    assertEquals(new HashingRegion(9, 1, 9, 34, 0, 34 + 33), getRegion(100, 33, RESULTS, 9));
    for (int i = 10; i < 100; i++) {
      assertEquals(HashingRegion.NONE, getRegion(100, 33, RESULTS, i));
    }
  }
  private static final MatchResult RESULTS_2 = initResults2();
  static MatchResult initResults2() {
    final MatchResult ret = new MatchResult(3);
    ret.addMatchResult(0, 1, 0, false);
    ret.addMatchResult(0, 50, 1, false);
    ret.addMatchResult(1, 32, 2, false);
    return ret;
  }

  public void testRegionsMultiMatch() {
    //3 threads
    assertEquals(new HashingRegion(0, 0, 0, 50, 0, 65), getRegion(3, 15, RESULTS_2, 0));
    AbstractAlignmentWriterThread.AlignmentWorkload work = getWorkload(3, 15, RESULTS_2, 0);
    assertEquals(0, work.getChunkStart());
    assertEquals(2, work.getChunkEnd());
    assertEquals(new HashingRegion(0, 50, 1, 32, 35, 47), getRegion(3, 15, RESULTS_2, 1));
    work = getWorkload(3, 15, RESULTS_2, 1);
    assertEquals(1, work.getChunkStart());
    assertEquals(3, work.getChunkEnd());
    assertEquals(new HashingRegion(1, 32, 1, 47, 17, 62), getRegion(3, 15, RESULTS_2, 2));
    work = getWorkload(3, 15, RESULTS_2, 2);
    assertEquals(2, work.getChunkStart());
    assertEquals(3, work.getChunkEnd());

    // 2 Threads
    assertEquals(new HashingRegion(0, 0, 0, 50, 0, 65), getRegion(2, 15, RESULTS_2, 0));
    assertEquals(new HashingRegion(0, 50, 1, 47, 35, 62), getRegion(2, 15, RESULTS_2, 1));
  }
  private static final MatchResult RESULTS_3 = initResult3();
  static MatchResult initResult3() {
    final MatchResult ret = new MatchResult(7);
    ret.addMatchResult(0, 1, 0, false);
    ret.addMatchResult(0, 10, 1, false);
    ret.addMatchResult(0, 50, 1, false);
    ret.addMatchResult(0, 60, 1, false);
    ret.addMatchResult(1, 32, 2, false);
    ret.addMatchResult(2, 42, 2, false);
    ret.addMatchResult(2, 52, 2, false);
    return ret;
  }

  public void testRegionsCloseMatches() {
    //3 threads
    assertEquals(new HashingRegion(0, 0, 0, 50, 0, 65), getRegion(3, 15, RESULTS_3, 0));
    AbstractAlignmentWriterThread.AlignmentWorkload work = getWorkload(3, 15, RESULTS_3, 0);
    assertEquals(0, work.getChunkStart());
    assertEquals(4, work.getChunkEnd());
    assertEquals(new HashingRegion(0, 50, 1, 32, 35, 47), getRegion(3, 15, RESULTS_3, 1));
    work = getWorkload(3, 15, RESULTS_3, 1);
    assertEquals(2, work.getChunkStart());
    assertEquals(5, work.getChunkEnd());
    assertEquals(new HashingRegion(1, 32, 2, 67, 17, 82), getRegion(3, 15, RESULTS_3, 2));
    work = getWorkload(3, 15, RESULTS_3, 2);
    assertEquals(4, work.getChunkStart());
    assertEquals(7, work.getChunkEnd());
  }

  static AbstractAlignmentWriterThread.AlignmentWorkload getWorkload(int nt, int padding, MatchResult results, int threadNo) {
    return new AbstractAlignmentWriterThread.AlignmentWorkload(nt, padding, results, threadNo);
  }

  static HashingRegion getRegion(int nt, int padding, MatchResult results, int threadNo) {
    final AbstractAlignmentWriterThread.AlignmentWorkload workload = new AbstractAlignmentWriterThread.AlignmentWorkload(nt, padding, results, threadNo);
    return workload.toRegion();
  }
}
