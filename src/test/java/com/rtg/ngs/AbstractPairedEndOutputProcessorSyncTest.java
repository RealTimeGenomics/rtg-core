/*
 * Copyright (c) 2018. Real Time Genomics Limited.
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
package com.rtg.ngs;

import java.io.File;
import java.io.IOException;
import java.util.HashSet;
import java.util.Set;

import com.rtg.index.hash.ngs.OutputProcessor;
import com.rtg.launcher.HashingRegion;
import com.rtg.launcher.SequenceParams;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.sam.RecordIterator;
import com.rtg.sam.SamRecordPopulator;
import com.rtg.sam.SamUtils;
import com.rtg.sam.ThreadedMultifileIterator;
import com.rtg.util.IORunnable;
import com.rtg.util.SimpleThreadPool;
import com.rtg.util.SingletonPopulatorFactory;
import com.rtg.util.StringUtils;
import com.rtg.util.TestUtils;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.TestDirectory;
import com.rtg.util.test.FileHelper;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.BlockCompressedInputStream;
import htsjdk.samtools.util.BlockCompressedInputStream.FileTermination;

/**
 *
 */
public abstract class AbstractPairedEndOutputProcessorSyncTest extends AbstractPairedEndOutputProcessorTest {

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
          mProc.process(0, "R", 1, 6, mThreadNum, 0);
          mProc.process(0, "F", 0, 2, mThreadNum, 0);
          mProc.process(0, "R", 1, 7, mThreadNum, 0);
          mProc.process(0, "F", 0, 3, mThreadNum, 0);
          mProc.process(0, "R", 1, 30, mThreadNum, 0);
          mProc.process(0, "F", 0, 20, mThreadNum, 0);
          mProc.process(0, "R", 1, 40, mThreadNum, 0);
        } finally {
          mProc.threadFinish();
        }
      } catch (final IOException e) {
        fail(e.getMessage());
      }
    }

  }

  public void testThreadClone1() throws Exception {
    checkThreadClone(null, false, false);
  }

  public void testThreadClone2() throws Exception {
    checkThreadClone(null, true, true);
  }

  public void testThreadClone3() throws Exception {
    try (TestDirectory tmp = new TestDirectory()) {
      checkThreadClone(tmp, false, true);
    }
  }

  public void testThreadCloneLocalZipText() throws Exception {
    try (TestDirectory tmp = new TestDirectory()) {
      checkThreadClone(tmp, true, false);
    }
  }

  abstract String getThreadCloneExpectedResource();

  abstract OutputFilter getOutputFilter();

  abstract OutputProcessor getPairedEndOutputProcessorSync(NgsParams param, MapStatistics stats, boolean outputUnmated, boolean outputUnmapped) throws IOException;

  abstract String getOutputFilePrefix();
  abstract String getOutputBamFilePrefix();

  /**
   * Test of threadClone method, of class PairedEndOutputProcessorSync.
   *
   * We use pairwise testing to systematically test the combinations of these
   * three parameters, each of which has two values.  That is, for each of the
   * three pairs of columns, we test all combinations of values for those columns.
   * (see <a href="http://pairwise.org">pairwise.org</a> for details).
   * <pre>
   *    checkThreadClone(null,     false, false);
   *    checkThreadClone(null,     true,  true);
   *    checkThreadClone(non-null, false, true);
   *    checkThreadClone(non-null, true,  false);
   * </pre>
   *
   * @param tempDir null or the name of a directory for temporary SAM files
   * @param gzipTemp true means compress the temporary SAM files
   * @param gzipResults true means compress the final result SAM files
   * @throws Exception on IO error
   */
  public void checkThreadClone(File tempDir, boolean gzipTemp, boolean gzipResults) throws Exception {
    final int numThreads = 4;
    final SimpleThreadPool stp = new SimpleThreadPool(numThreads, "TestPairedEnd", true);
    final NgsParams params = getDefaultBuilder(tempDir, gzipResults, getOutputFilter(), null).numberThreads(numThreads).unknownsPenalty(9).create();
    try (OutputProcessor sync = getPairedEndOutputProcessorSync(params, null, true, true)) {
      for (int i = 0; i < numThreads; ++i) {
        final long padding = params.calculateThreadPadding();
        final long start = i * MAX_COORD / numThreads;
        final long end = (i + 1) * MAX_COORD / numThreads;
        final HashingRegion region = new HashingRegion(0, start, 0, end, Math.max(0, start - padding), Math.min(MAX_COORD, end + padding));
        //System.err.println(region);
        stp.execute(new SimpleProcess(sync, region, i));
      }
      stp.terminate();
      sync.finish();
    } finally {
      params.close();
    }
    final String contents = gzipResults
        ? TestUtils.stripSAMHeader(FileHelper.gzFileToString(new File(new File(mDir, "hitDir"), getOutputFilePrefix() + FileUtils.GZ_SUFFIX)))
            : TestUtils.stripSAMHeader(FileUtils.fileToString(new File(new File(mDir, "hitDir"), getOutputFilePrefix())));
        mNano.check(getThreadCloneExpectedResource(), contents, false);
  }


  static class SimpleProcess2 implements IORunnable {

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

  static final String SORTING_TEMPLATE = ">t" + StringUtils.LS + "tgcaagacaagagggcctcc"
      + "tgcaagacaagagggcctcc"
      + "tgcaagacaagagggcctcc" + StringUtils.LS;

  static final String SORTING_READ_LEFT = ">r" + StringUtils.LS + TEMP_LEFT + StringUtils.LS
      + ">s" + StringUtils.LS + TEMP_LEFT + StringUtils.LS
      + ">t" + StringUtils.LS + TEMP_LEFT + StringUtils.LS;

  static final String SORTING_READ_RIGHT = ">r" + StringUtils.LS + TEMP_RIGHT + StringUtils.LS
      + ">s" + StringUtils.LS + TEMP_RIGHT + StringUtils.LS
      + ">t" + StringUtils.LS + TEMP_RIGHT + StringUtils.LS;

  private NgsParamsBuilder getSortedAppendBuilder(File tempDir, boolean gzipOutputs, boolean bam) throws IOException {
    final File templateok = FileUtils.createTempDir("template", "ngs", mDir);
    final File leftok = FileUtils.createTempDir("left", "ngs", mDir);
    final File rightok = FileUtils.createTempDir("right", "ngs", mDir);
    final File hitsDir = new File(mDir, "hitDir");

    ReaderTestUtils.getReaderDNA(SORTING_TEMPLATE, templateok, null).close();
    ReaderTestUtils.getReaderDNA(SORTING_READ_LEFT, leftok, null).close();
    ReaderTestUtils.getReaderDNA(SORTING_READ_RIGHT, rightok, null).close();

    final NgsFilterParams filterParams = NgsFilterParams.builder().outputFilter(getOutputFilter())
        .topN(10).errorLimit(5).zip(gzipOutputs).create();
    final NgsOutputParams outputParams = NgsOutputParams.builder()
        .progress(false).outputDir(hitsDir)
        .tempFilesDir(tempDir).filterParams(filterParams).bam(bam).create();

    return NgsParams.builder()
        .buildFirstParams(SequenceParams.builder().directory(leftok).useMemReader(true).create())
        .buildSecondParams(SequenceParams.builder().directory(rightok).useMemReader(true).create())
        .searchParams(SequenceParams.builder().directory(templateok).useMemReader(true).loadNames(true).create())
        .outputParams(outputParams)
        .maskParams(new NgsMaskParamsGeneral(4, 1, 1, 1))
        .substitutionPenalty(1).gapOpenPenalty(1).gapExtendPenalty(1).unknownsPenalty(0)
        .maxFragmentLength(1000).minFragmentLength(0);
  }

  void sortedAppendTest(String expectedResource) throws IOException {
    final File dir = FileUtils.createTempDir("dop", "test");
    try {
      final NgsParams params = getSortedAppendBuilder(dir, true, false).numberThreads(4).create();
      final OutputProcessor sync = getPairedEndOutputProcessorSync(params, null, true, true);
      final OutputProcessor one = sync.threadClone(new HashingRegion(0, 20, 0, 40, -1, -1));
      final OutputProcessor two = sync.threadClone(new HashingRegion(0, 40, 0, 60, -1, -1));
      final OutputProcessor three = sync.threadClone(new HashingRegion(0, 0, 0, 20, -1, -1));
      one.process(0, "F", 2, 20, 0, 0);
      one.process(0, "R", 3, 25, 0, 0);
      two.process(0, "F", 4, 40, 0, 0);
      two.process(0, "R", 5, 45, 0, 0);
      three.process(0, "F", 0, 0, 0, 0);
      three.process(0, "R", 1, 5, 0, 0);
      one.threadFinish();
      two.threadFinish();
      three.threadFinish();
      sync.finish();
      final File gzipFile = new File(new File(mDir, "hitDir"), getOutputFilePrefix() + FileUtils.GZ_SUFFIX);
      assertTrue(gzipFile.exists());
      assertEquals(FileTermination.HAS_TERMINATOR_BLOCK, BlockCompressedInputStream.checkTermination(gzipFile));
      final String contents = TestUtils.stripSAMHeader(FileHelper.gzFileToString(gzipFile));
      mNano.check(expectedResource, contents, false);
    } finally {
      assertTrue(FileHelper.deleteAll(dir));
    }
  }

  void sortedAppendTestBam(String expectedResource) throws IOException {
    final File dir = FileUtils.createTempDir("dop", "test");
    try {
      final NgsParams params = getSortedAppendBuilder(dir, true, true).numberThreads(4).create();
      final OutputProcessor sync = getPairedEndOutputProcessorSync(params, null, true, true);
      final OutputProcessor one = sync.threadClone(new HashingRegion(0, 20, 0, 40, -1, -1));
      final OutputProcessor two = sync.threadClone(new HashingRegion(0, 40, 0, 60, -1, -1));
      final OutputProcessor three = sync.threadClone(new HashingRegion(0, 0, 0, 20, -1, -1));
      one.process(0, "F", 2, 20, 0, 0);
      one.process(0, "R", 3, 25, 0, 0);
      two.process(0, "F", 4, 40, 0, 0);
      two.process(0, "R", 5, 45, 0, 0);
      three.process(0, "F", 0, 0, 0, 0);
      three.process(0, "R", 1, 5, 0, 0);
      one.threadFinish();
      two.threadFinish();
      three.threadFinish();
      sync.finish();
      final File bamFile = new File(new File(mDir, "hitDir"), getOutputBamFilePrefix());
      assertTrue(bamFile.exists());
      assertTrue(SamUtils.isBAMFile(bamFile));
      assertEquals(FileTermination.HAS_TERMINATOR_BLOCK, BlockCompressedInputStream.checkTermination(bamFile));

      final Set<File> files = new HashSet<>();
      files.add(bamFile);
      final StringBuilder sb = new StringBuilder();
      try (RecordIterator<SAMRecord> mfi = new ThreadedMultifileIterator<>(files, new SingletonPopulatorFactory<>(new SamRecordPopulator()))) {
        while (mfi.hasNext()) {
          final SAMRecord rec = mfi.next();
          sb.append(rec.getSAMString());
        }
      } finally {
        final String contents = sb.toString();
        mNano.check(expectedResource, contents, false);
      }
    } finally {
      assertTrue(FileHelper.deleteAll(dir));
    }
  }

}
