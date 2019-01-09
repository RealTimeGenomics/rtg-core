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

import com.rtg.bed.BedUtils;
import com.rtg.calibrate.CalibratorTest;
import com.rtg.index.hash.ngs.OutputProcessor;
import com.rtg.launcher.CommonFlags;
import com.rtg.launcher.DefaultReaderParams;
import com.rtg.launcher.HashingRegion;
import com.rtg.launcher.SequenceParams;
import com.rtg.mode.SequenceType;
import com.rtg.ngs.tempstage.PairedTempFileWriterImpl;
import com.rtg.pairedend.SlidingWindowCollector;
import com.rtg.reader.MockSequencesReader;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.reader.SequencesReader;
import com.rtg.tabix.TabixIndexer;
import com.rtg.util.InvalidParamsException;
import com.rtg.util.SimpleThreadPool;
import com.rtg.util.StringUtils;
import com.rtg.util.TestUtils;
import com.rtg.util.cli.CommandLine;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.util.io.TestDirectory;
import com.rtg.util.test.FileHelper;

import htsjdk.samtools.SAMReadGroupRecord;

/**
 */
public class TopNPairedEndOutputProcessorSyncTest extends AbstractPairedEndOutputProcessorSyncTest {

  @Override
  OutputProcessor getPairedEndOutputProcessorSync(NgsParams param, MapStatistics stats, boolean outputUnmated, boolean outputUnmapped) throws IOException {
    return new TopNPairedEndOutputProcessorSync(param, stats, outputUnmated, outputUnmapped);
  }

  @Override
  PairedEndOutputProcessor getPEOP(PairedTempFileWriterImpl writer, SlidingWindowCollector collector) {
    return new PairedEndOutputProcessor(writer, collector);
  }

  @Override
  String getThreadCloneExpectedResource() {
    return "tnpeops_tc.txt";
  }

  @Override
  OutputFilter getOutputFilter() {
    return OutputFilter.TOPN_PAIRED_END;
  }

  @Override
  String getOutputFilePrefix() {
    return NgsOutputParams.MATED_SAM_FILE_NAME;
  }
  @Override
  String getOutputBamFilePrefix() {
    return NgsOutputParams.MATED_BAM_FILE_NAME;
  }

  public void testUnmatedLocal() throws IOException, InvalidParamsException {
    try (TestDirectory tmp = new TestDirectory()) {
      checkUnmated(tmp);
    }
  }

  public void testUnmatedRemote() throws IOException, InvalidParamsException {
    checkUnmated(null);
  }

  /**
   * Test of threadClone method, of class PairedEndOutputProcessorSync.
   */
  public void checkUnmated(final File tempDir) throws IOException, InvalidParamsException {
    final ByteArrayOutputStream logBytes = new ByteArrayOutputStream();
    try (PrintStream log = new PrintStream(logBytes)) {
      Diagnostic.setLogStream(log);
      final File rgFile = FileHelper.resourceToFile("com/rtg/sam/resources/readgroup_cg.txt", new File(tempDir, "readGroup.txt"));
      try {
        final int numThreads = 4;
        final SimpleThreadPool stp = new SimpleThreadPool(numThreads, "TestUnmated", true);
        final NgsParams params = getDefaultBuilder(tempDir, false, OutputFilter.TOPN_PAIRED_END, rgFile).numberThreads(numThreads).create();
        final SequencesReader ref = params.searchParams().reader();
        final HashingRegion[] regions = HashingRegion.splitWorkload(ref, params.sex(), 0, ref.numberSequences(), params.numberThreads() * params.threadMultiplier(), HashingRegion.DEFAULT_MIN_CHUNK_SIZE, params.calculateThreadPadding());
        try (TopNPairedEndOutputProcessorSync sync = new TopNPairedEndOutputProcessorSync(params, null, true, true)) {
          for (int i = 0; i < regions.length; ++i) {
            stp.execute(new SimpleProcess2(sync, regions[i], i));
          }
          stp.terminate();
          sync.finish();
        }
        final String contentsUnmapped = TestUtils.stripSAMHeader(FileUtils.fileToString(new File(new File(mDir, "hitDir"), NgsOutputParams.UNMAPPED_SAM_FILE_NAME)));
        final String contentsUnmated = TestUtils.stripSAMHeader(FileUtils.fileToString(new File(new File(mDir, "hitDir"), NgsOutputParams.UNMATED_SAM_FILE_NAME)));
        final String contentsUnmatedCal = CalibratorTest.stripVersion(FileUtils.fileToString(new File(new File(mDir, "hitDir"), NgsOutputParams.UNMATED_SAM_FILE_NAME + CommonFlags.RECALIBRATE_EXTENSION)));
        mNano.check("topnpeops-checkunmated-unmated", contentsUnmated, false);
        mNano.check("topnpeops-checkunmated-unmapped", contentsUnmapped, false);
        mNano.check("topnpeops-checkunmated-unmated-cal", contentsUnmatedCal, true);
      } finally {
        assertTrue(FileHelper.deleteAll(rgFile));
      }
    }
  }

  static final String TEMPLATE_LARGE = ">t" + StringUtils.LS + "TGCAAGACAAGAGGGCCTCCTTCTTTCAGTGATGCTGATCGACTCACCAGCTAGCTCAGCATACTCAGACCACTTGCTAGGATCGGCTAAAATCGCTGTCTCGATCGATAATCGATGGCTCGGGATCGA" + StringUtils.LS;
  static final String READS_LEFT = ">r1" + StringUtils.LS + "TGCAAGACAAGAGGGCCTCC" + StringUtils.LS
                                 + ">r2" + StringUtils.LS + "GTGATGCTGATCGACTCACC" + StringUtils.LS
                                 + ">r3" + StringUtils.LS + "ACAAGAGGGCCTCCTTCTTT" + StringUtils.LS;
  static final String READS_RIGHT = ">r1" + StringUtils.LS + "TCACCAGCTAGCTCAGCATA" + StringUtils.LS
                                  + ">r2" + StringUtils.LS + "CAGATGTATTGGATTTAGGG" + StringUtils.LS
                                  + ">r3" + StringUtils.LS + "AATCGATGGCTCGGGATCGA" + StringUtils.LS;

  public void testAugmented() throws IOException, InvalidParamsException {
    final ByteArrayOutputStream logBytes = new ByteArrayOutputStream();
    try (TestDirectory tmp = new TestDirectory(); PrintStream log = new PrintStream(logBytes)) {
      Diagnostic.setLogStream(log);
      final File templateok = new File(tmp, "template");
      final File input = new File(tmp, "reads");
      assertTrue(input.mkdir());
      final File leftok = new File(input, "left");
      final File rightok = new File(input, "right");
      ReaderTestUtils.getReaderDNA(TEMPLATE_LARGE, templateok, null).close();
      ReaderTestUtils.getReaderDNA(READS_LEFT, leftok, null).close();
      ReaderTestUtils.getReaderDNA(READS_RIGHT, rightok, null).close();
      final File outDir = new File(tmp, "outDir");
      new MapCli().mainInit(new String[]{"-t", templateok.getPath(), "-i", input.getPath(), "-o", outDir.getPath(), "--sam-rg", "@RG\\tID:RG23\\tSM:NA123\\tPL:COMPLETE", "-Z", "--max-fragment-size", "65", "--sam", "--no-merge"}, logBytes, log);
      final String contentsUnmapped = TestUtils.stripSAMHeader(FileUtils.fileToString(new File(outDir, NgsOutputParams.UNMAPPED_SAM_FILE_NAME)));
      final String contentsUnmated = TestUtils.stripSAMHeader(FileUtils.fileToString(new File(outDir, NgsOutputParams.UNMATED_SAM_FILE_NAME)));

      mNano.check("topnpeops-checkunmatedaugment-unmated", contentsUnmated, false);
      mNano.check("topnpeops-checkunmatedaugment-unmapped", contentsUnmapped, false);
    }
  }

  public void testSortedAppend() throws IOException {
    sortedAppendTest("tnpeops_sortedAppend.txt");
  }
  public void testSortedAppendBam() throws IOException {
    sortedAppendTestBam("tnpeops_sortedAppend.txt");
  }

  public void test() throws IOException {
    try (TestDirectory dir = new TestDirectory()) {
      final File foo = new File(dir, "foo");
      final File fooTbi = new File(dir, "foo.tbi");
      final File fooBai = new File(dir, "foo.bai");
      assertEquals(fooTbi, AbstractMulticoreFilterConcat.indexFileName(foo, false));
      assertEquals(fooBai, AbstractMulticoreFilterConcat.indexFileName(foo, true));
    }
  }

  private static final String BED = ""
      + "t\t50\t200\tfail"
      ;
  public void testCalibrateRegion() throws IOException, InvalidParamsException {
    try (TestDirectory tempDir = new TestDirectory()) {
      final File bedFile = new File(tempDir, "bed");
      FileUtils.stringToFile(BED, bedFile);
      CommandLine.clearCommandArgs();
      final ByteArrayOutputStream logBytes = new ByteArrayOutputStream();
      try (PrintStream log = new PrintStream(logBytes)) {
        Diagnostic.setLogStream(log);
        final File rgFile = FileHelper.resourceToFile("com/rtg/sam/resources/readgroup_cg.txt", new File(tempDir, "readGroup.txt"));
        try {
          final int numThreads = 4;
          final SimpleThreadPool stp = new SimpleThreadPool(numThreads, "TestUnmated", true);
          final File hitDir = new File(mDir, "hitDir");
          final NgsOutputParams output = NgsOutputParams.builder().calibrate(true).calibrateRegions(BedUtils.regions(bedFile)).readGroup(new SAMReadGroupRecord("RG23")).outputDir(hitDir).create();
          final NgsParams params = getDefaultBuilder(tempDir, false, OutputFilter.TOPN_PAIRED_END, rgFile).numberThreads(numThreads).outputParams(output).create();
          final SequencesReader ref = params.searchParams().reader();
          final HashingRegion[] regions = HashingRegion.splitWorkload(ref, params.sex(), 0, ref.numberSequences(), params.numberThreads() * params.threadMultiplier(), HashingRegion.DEFAULT_MIN_CHUNK_SIZE, params.calculateThreadPadding());
          try (TopNPairedEndOutputProcessorSync sync = new TopNPairedEndOutputProcessorSync(params, null, true, true)) {
            for (int i = 0; i < regions.length; ++i) {
              stp.execute(new SimpleProcess2(sync, regions[i], i));
            }
            stp.terminate();
            sync.finish();
          }
          final String contentsUnmapped = TestUtils.stripSAMHeader(FileUtils.fileToString(new File(hitDir, NgsOutputParams.UNMAPPED_SAM_FILE_NAME)));
          final String contentsUnmated = TestUtils.stripSAMHeader(FileUtils.fileToString(new File(hitDir, NgsOutputParams.UNMATED_SAM_FILE_NAME)));
          final String calibrateFile = FileUtils.fileToString(new File(hitDir, NgsOutputParams.UNMATED_SAM_FILE_NAME + CommonFlags.RECALIBRATE_EXTENSION));
          final String contentsUnmatedCal = CalibratorTest.stripVersion(calibrateFile);
          mNano.check("topnpeops-checkunmated-unmated", contentsUnmated, false);
          mNano.check("topnpeops-checkunmated-unmapped", contentsUnmapped, false);
          mNano.check("topnpeops-checkunmated-unmated-region-cal", contentsUnmatedCal, true);
        } finally {
          assertTrue(FileHelper.deleteAll(rgFile));
        }
      }
    }
  }

  public void testLargeSequence() throws IOException {
    try (TestDirectory dir = new TestDirectory()) {
      final MockSequencesReader r = new MockSequencesReader(SequenceType.DNA, 10) {
        @Override
        public long maxLength() {
          return TabixIndexer.MAXIMUM_REFERENCE_LENGTH + 1;
        }
      };
      final MemoryPrintStream mps = new MemoryPrintStream();
      Diagnostic.setLogStream(mps.printStream());
      final SequenceParams sequenceParams = new SequenceParams.SequenceParamsBuilder().readerParam(new DefaultReaderParams(r)).create();
      final NgsOutputParams outputParams = new NgsOutputParamsBuilder().outputIndex(true).bam(false).create();
      final NgsParams params = new NgsParamsBuilder().searchParams(sequenceParams).outputParams(outputParams).create();
      final SansFilterConcat sans = new SansFilterConcat(params, new ReadStatusTracker(10, new PairedEndMapStatistics(false, null)), 0);
      final AbstractMulticoreFilterConcat.OutputWrapper streams = sans.createStreams(1, new File[]{new File(dir, "Foo")}, new File[]{new File(dir, "bar")}, false, true, 0);
      try {
        final String s = mps.toString();
        assertTrue(s.contains("Cannot produce TABIX index"));
        assertTrue(s.contains("as maximum reference sequence length is exceeded."));
      } finally {
        streams.mOutputStream.close();
      }
    }
  }

}
