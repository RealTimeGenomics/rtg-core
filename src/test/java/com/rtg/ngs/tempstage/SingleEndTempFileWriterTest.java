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
package com.rtg.ngs.tempstage;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.PrintStream;

import com.rtg.launcher.HashingRegion;
import com.rtg.launcher.SequenceParams;
import com.rtg.ngs.MapStatistics;
import com.rtg.reader.Arm;
import com.rtg.ngs.MapStatisticsField;
import com.rtg.ngs.NgsFilterParams;
import com.rtg.ngs.NgsOutputParams;
import com.rtg.ngs.NgsParams;
import com.rtg.ngs.OutputFilter;
import com.rtg.ngs.ReadStatusTracker;
import com.rtg.ngs.SharedResources;
import com.rtg.ngs.SingleEndMapStatistics;
import com.rtg.ngs.blocking.MapQScoringReadBlocker;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.util.IntegerOrPercentage;
import com.rtg.util.MaxShiftFactor;
import com.rtg.util.StringUtils;
import com.rtg.util.TestUtils;
import com.rtg.util.cli.CommandLine;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.FileUtils;
import com.rtg.util.test.FileHelper;

import htsjdk.samtools.SAMReadGroupRecord;

import junit.framework.TestCase;

/**
 * Tests corresponding class.
 */
public class SingleEndTempFileWriterTest extends TestCase {

  private File mDir;

  static final String LS = StringUtils.LS;

  @Override
  public void setUp() throws Exception {
    mDir = FileHelper.createTempDirectory();
    Diagnostic.setLogStream();
  }

  @Override
  public void tearDown() {
    FileHelper.deleteAll(mDir);
    mDir = null;
  }

  static final String TEMPLATE = ">t" + StringUtils.LS + "acacactgcaagacaagagggcctcccacagcactctcagcccacactggtcgggggccaaagggg" + StringUtils.LS;
  static final String TEMP1 = "acacactgcaagcaagagggcctccc";
  static final String READS = ">r" + StringUtils.LS + TEMP1 + StringUtils.LS;

  public void testAlignment() throws IOException {
    final File template = FileUtils.createTempDir("template", "ngs", mDir);
    final File inDir = FileUtils.createTempDir("reads", "ngs", mDir);

    ReaderTestUtils.getReaderDNA(TEMPLATE, template, null).close();
    ReaderTestUtils.getReaderDNA(READS, inDir, null).close();

    CommandLine.setCommandArgs("wibble", "-h", "dream-turnip");
    try {

      final SAMReadGroupRecord rg = new SAMReadGroupRecord("L10");
      rg.setSample("NA1234");
      final NgsOutputParams op = NgsOutputParams.builder()
              .filterParams(NgsFilterParams.builder().outputFilter(OutputFilter.SAM_SINGLE_END).create())
              .readGroup(rg)
              .create();
      final NgsParams params = NgsParams.builder()
              .buildFirstParams(SequenceParams.builder().directory(inDir).useMemReader(true).create())
              .searchParams(SequenceParams.builder().directory(template).useMemReader(true).loadNames(true).create())
              .maxFragmentLength(1000).minFragmentLength(0)
              .outputParams(op)
              .gapOpenPenalty(1).gapExtendPenalty(1).substitutionPenalty(1).unknownsPenalty(0)
              .create();
      final ByteArrayOutputStream log = new ByteArrayOutputStream();
      final PrintStream prLog = new PrintStream(log);
      Diagnostic.setLogStream(prLog);


      final MapQScoringReadBlocker blocker = new MapQScoringReadBlocker((int) params.buildFirstParams().reader().numberSequences(), 10);

      try (ByteArrayOutputStream outStream = new ByteArrayOutputStream()) {
        try (SingleEndTempFileWriter w = new SingleEndTempFileWriter(params, new PairedTempFileWriterImplTest.UselessStatusIdListener(), SharedResources.generateSharedResources(params))) {
          w.initialiseAlignments(outStream, blocker);
          w.nextTemplateId(0);
          final boolean result = w.alignmentResult(0, false, 0);
          assertTrue(result);
          final BinaryTempFileRecord record = w.alignmentResultUnfiltered(0, false, 0);
          assertEquals(1, record.getStartPosition());
          assertEquals("12=1D14=", new String(record.getCigarString()));
          //assertEquals("L10", record.getAttribute("RG")); no longer included in temporary file
          //assertEquals(1, record.getHeader().getReadGroups().size());
          //assertEquals("L10", record.getHeader().getReadGroups().get(0).getReadGroupId());
          //assertEquals("NA1234", record.getHeader().getReadGroups().get(0).getSample());
          assertNotNull(w.getBlocker());
        }
      } finally {
        Diagnostic.setLogStream();
        prLog.close();

      }
      TestUtils.containsAll(log.toString(),
              "Statistics of multithreaded blocked",
              "1 reads had count 0",
              "Total reads 1",
              "Alignment score filter ",
              ": 2 passed, 2 total",
              "Duplicates detected during SAM writing: 1");
      assertFalse(log.toString().contains("Reordering buffer overflow caused "));
    } finally {
      CommandLine.clearCommandArgs();
    }
  }


  private class MySamSingleEndAlignmentWriter extends SingleEndTempFileWriter {

    final int[] mMyMatchActions;
    MySamSingleEndAlignmentWriter(NgsParams param, ReadStatusTracker stat, int[] myMatchActions)
    throws IOException {
      super(param, stat, SharedResources.generateSharedResources(param));
      mMyMatchActions = myMatchActions;
    }

    @Override
    protected int[] calculateEditDistance(byte[] read, int length, int start, boolean rc, IntegerOrPercentage maxMismatches, boolean cgLeft, int readId) {
      return mMyMatchActions;
    }
  }

  private class MyStatusReadIdListenerImpl extends ReadStatusTracker {
    MyStatusReadIdListenerImpl(int numReads, MapStatistics stats) {
      super(numReads, stats);
    }

    @Override
    protected void calculateStatistics(boolean pairedEnd, boolean allhits) {
      super.calculateStatistics(pairedEnd, allhits);
    }
  }

  static final String TSTR = "acgt";
  static final String READS2 = ">r1" + LS + TSTR + TSTR + TSTR + TSTR + TSTR + TSTR
  + ">r2" + LS + TSTR + TSTR + TSTR + TSTR + TSTR + TSTR
  + ">r3" + LS + TSTR + TSTR + TSTR + TSTR + TSTR + TSTR
  + ">r4" + LS + TSTR + TSTR + TSTR + TSTR + TSTR + TSTR
  + ">r5" + LS + TSTR + TSTR + TSTR + TSTR + TSTR + TSTR
  + LS;

  static final int[] ACTIONS1 = {99, 10, 5, 2};


  public void testAlignment2() throws IOException {
    final File template = FileUtils.createTempDir("template", "ngs", mDir);
    final File inDir = FileUtils.createTempDir("reads", "ngs", mDir);

    ReaderTestUtils.getReaderDNA(TEMPLATE, template, null).close();
    ReaderTestUtils.getReaderDNA(READS2, inDir, null).close();

    final NgsOutputParams op = NgsOutputParams.builder()
    .filterParams(NgsFilterParams.builder().outputFilter(OutputFilter.SAM_SINGLE_END).create()).create();
    final NgsParams params = NgsParams.builder()
    .buildFirstParams(SequenceParams.builder().directory(inDir).useMemReader(true).create())
    .searchParams(SequenceParams.builder().directory(template).useMemReader(true).loadNames(true).create())
    .maxFragmentLength(1000).minFragmentLength(0)
    .outputParams(op)
    .create();
    final ByteArrayOutputStream log = new ByteArrayOutputStream();
    final PrintStream prLog = new PrintStream(log);
    Diagnostic.setLogStream(prLog);


    final MapQScoringReadBlocker blocker = new MapQScoringReadBlocker((int) params.buildFirstParams().reader().numberSequences(), 10);

    try (ByteArrayOutputStream outStream = new ByteArrayOutputStream()) {
      final MapStatistics stats = new SingleEndMapStatistics(null);
      final MyStatusReadIdListenerImpl listener = new MyStatusReadIdListenerImpl(1, stats);
      try (MySamSingleEndAlignmentWriter w = new MySamSingleEndAlignmentWriter(params, listener, ACTIONS1)) {
        final HashingRegion region = new HashingRegion(0, 1, 0, 99, -1, -1);
        w.setClipRegion(region);
        w.initialiseAlignments(outStream, blocker);
        w.nextTemplateId(0);
        boolean result = w.alignmentResult(0, false, 1);
        assertTrue(result);
        result = w.alignmentResult(0, false, 11);
        assertTrue(result);
        w.alignmentResultUnfiltered(0, false, 0);
        listener.calculateStatistics(false, false);
        assertNotNull(w.getBlocker());
      }
    } finally {
      Diagnostic.setLogStream();
      prLog.close();

    }
    TestUtils.containsAll(log.toString(),
        "Statistics of multithreaded blocked",
        "5 reads had count 0",
        "Total reads 5",
        "Alignment score filter ",
        ": 0 passed, 0 total",
        "Duplicates detected during SAM writing: 0");
  }

  public void testAlignment3() throws IOException {
    final File template = FileUtils.createTempDir("template", "ngs", mDir);
    final File inDir = FileUtils.createTempDir("reads", "ngs", mDir);

    ReaderTestUtils.getReaderDNA(TEMPLATE, template, null).close();
    ReaderTestUtils.getReaderDNA(READS2, inDir, null).close();

    final NgsOutputParams op = NgsOutputParams.builder()
      .filterParams(NgsFilterParams.builder().outputFilter(OutputFilter.SAM_SINGLE_END).matedMaxMismatches(IntegerOrPercentage.valueOf(15)).create()).create();
    final NgsParams params = NgsParams.builder()
    .buildFirstParams(SequenceParams.builder().directory(inDir).useMemReader(true).create())
    .searchParams(SequenceParams.builder().directory(template).useMemReader(true).loadNames(true).create())
    .maxFragmentLength(1000).minFragmentLength(0)
    .outputParams(op)
    .gapOpenPenalty(1).gapExtendPenalty(1).substitutionPenalty(1).unknownsPenalty(1)
    .alignerBandWidthFactor(new MaxShiftFactor(1))
    .create();
    final ByteArrayOutputStream log = new ByteArrayOutputStream();
    final PrintStream prLog = new PrintStream(log);
    Diagnostic.setLogStream(prLog);

    //    final ByteArrayOutputStream outUnmappedStream = new ByteArrayOutputStream();
    final MapQScoringReadBlocker blocker = new MapQScoringReadBlocker((int) params.buildFirstParams().reader().numberSequences(), 10);
    try (ByteArrayOutputStream outStream = new ByteArrayOutputStream()) {
      final MapStatistics stats = new SingleEndMapStatistics(null);
      final MyStatusReadIdListenerImpl listener = new MyStatusReadIdListenerImpl(5, stats);
      try (SingleEndTempFileWriter w = new SingleEndTempFileWriter(params, listener, SharedResources.generateSharedResources(params))) {
        final HashingRegion region = new HashingRegion(0, 1, 0, 100, -1, -1);
        w.setClipRegion(region);

        w.initialiseAlignments(outStream, blocker);
//        w.initialiseUnmapped(outUnmappedStream, false);
        w.nextTemplateId(0);
        boolean result;
        w.setClipRegion(new HashingRegion(1, 2));
        result = w.alignmentResult(0, false, 1000);
        for (int i = 0; i < 10; ++i) {
          result = w.alignmentResult(0, false, 1 + 4 * i);
        }
        assertTrue(result);
        w.setClipRegion(new HashingRegion(0, 1, 0, 999, -1, -1));
        result = w.alignmentResult(1, false, 1000);
        assertTrue(result);
        w.setClipRegion(new HashingRegion(0, 2));
        for (int i = 0; i < 10; ++i) {
          result = w.alignmentResult(0, false, 1 + 4 * i);
        }
        assertFalse(result);

        final BinaryTempFileRecord record = w.alignmentResultUnfiltered(0, false, 0);
        assertEquals(1, record.getStartPosition());
        assertEquals("2=2X2=6X2=4X1=2X1=1X1=", new String(record.getCigarString()));
        listener.calculateStatistics(false, false);
//        TestUtils.containsAll(StringUtils.toString(stats), "UNMAPPED_READS=4", "MISSING_LEFT=0", "UNMAPPED_READS_PERCENT=80",
//            "MATED_UNIQUE_READS=0", "MATED_AMBIG_READS=0", "MISSING_RIGHT=0", "UNMATED_AMBIG_READS_PERCENT=20", "UNMATED_AMBIG_READS=1", "MATED_UNIQUE_READS_PERCENT=0",
//            "TOTAL_READS=5", "MAPPED_READS_PERCENT=20", "MAPPED_READS=1");
        assertEquals(4L, stats.value(MapStatisticsField.UNMAPPED_NO_HITS, Arm.LEFT));
        assertEquals(0L, stats.value(MapStatisticsField.MISSING, Arm.LEFT));
        assertEquals(0L, stats.value(MapStatisticsField.MATED_UNIQUE_READS, Arm.LEFT));
        assertEquals(0L, stats.value(MapStatisticsField.MATED_AMBIG_READS, Arm.LEFT));
        assertEquals(1L, stats.value(MapStatisticsField.UNMATED_AMBIG_READS, Arm.LEFT));
        assertEquals(5L, stats.value(MapStatisticsField.TOTAL_READS, Arm.LEFT));
        assertEquals(80.0, stats.valueAsPercent(MapStatisticsField.UNMAPPED_NO_HITS, Arm.LEFT));
        assertEquals(0.0, stats.valueAsPercent(MapStatisticsField.MISSING, Arm.LEFT));
        assertEquals(0.0, stats.valueAsPercent(MapStatisticsField.MATED_UNIQUE_READS, Arm.LEFT));
        assertEquals(0.0, stats.valueAsPercent(MapStatisticsField.MATED_AMBIG_READS, Arm.LEFT));
        assertEquals(20.0, stats.valueAsPercent(MapStatisticsField.UNMATED_AMBIG_READS, Arm.LEFT));
        assertEquals(100.0, stats.valueAsPercent(MapStatisticsField.TOTAL_READS, Arm.LEFT));

        stats.reset();

        w.alignmentResult(1, true, 11);
        w.alignmentResult(2, true, 26);
//        w.unmappedResult(3, (char) 0);
        w.alignmentResult(4, true, 20);
        w.alignmentResult(4, true, 120);
        listener.calculateStatistics(false, false);
//        TestUtils.containsAll(StringUtils.toString(stats), "UNMAPPED_READS=4", "MISSING_LEFT=0", "UNMAPPED_READS_PERCENT=80",
//           "MATED_UNIQUE_READS=0", "MATED_AMBIG_READS=0", "MISSING_RIGHT=0", "UNMATED_AMBIG_READS_PERCENT=20", "UNMATED_AMBIG_READS=1", "MATED_UNIQUE_READS_PERCENT=0",
//            "TOTAL_READS=5", "MAPPED_READS_PERCENT=20", "MAPPED_READS=1");
        assertEquals(4L, stats.value(MapStatisticsField.UNMAPPED_NO_HITS, Arm.LEFT));
        assertEquals(0L, stats.value(MapStatisticsField.MISSING, Arm.LEFT));
        assertEquals(0L, stats.value(MapStatisticsField.MATED_UNIQUE_READS, Arm.LEFT));
        assertEquals(0L, stats.value(MapStatisticsField.MATED_AMBIG_READS, Arm.LEFT));
        assertEquals(1L, stats.value(MapStatisticsField.UNMATED_AMBIG_READS, Arm.LEFT));
        assertEquals(5L, stats.value(MapStatisticsField.TOTAL_READS, Arm.LEFT));
        assertEquals(80.0, stats.valueAsPercent(MapStatisticsField.UNMAPPED_NO_HITS, Arm.LEFT));
        assertEquals(0.0, stats.valueAsPercent(MapStatisticsField.MISSING, Arm.LEFT));
        assertEquals(0.0, stats.valueAsPercent(MapStatisticsField.MATED_UNIQUE_READS, Arm.LEFT));
        assertEquals(0.0, stats.valueAsPercent(MapStatisticsField.MATED_AMBIG_READS, Arm.LEFT));
        assertEquals(20.0, stats.valueAsPercent(MapStatisticsField.UNMATED_AMBIG_READS, Arm.LEFT));
        assertEquals(100.0, stats.valueAsPercent(MapStatisticsField.TOTAL_READS, Arm.LEFT));
        assertNotNull(w.getBlocker());
      }
    } finally {
      Diagnostic.setLogStream();
      prLog.close();

    }
    TestUtils.containsAll(log.toString(),
        "Statistics of multithreaded blocked",
        "5 reads had count 0",
        "Total reads 5",
        "Alignment score filter ",
        "Duplicates detected during SAM writing: 6");
  }

  public void testAlignment4() throws IOException {
    final File template = FileUtils.createTempDir("template", "ngs", mDir);
    final File inDir = FileUtils.createTempDir("reads", "ngs", mDir);

    final StringBuilder reads3 = new StringBuilder();
    for (int i = 1; i < 2051; ++i) {
      reads3.append(">r").append(i).append(LS).append("cgt").append(TSTR).append(TSTR).append(TSTR).append(TSTR).append(TSTR);
    }
    for (int i = 2051; i < 2100; ++i) {
      reads3.append(">r").append(i).append(LS).append(TSTR).append(TSTR).append(TSTR).append(TSTR).append(TSTR).append(TSTR);
    }

    ReaderTestUtils.getReaderDNA(TEMPLATE, template, null).close();
    ReaderTestUtils.getReaderDNA(reads3.toString(), inDir, null).close();

    final NgsOutputParams op = NgsOutputParams.builder()
      .filterParams(NgsFilterParams.builder().outputFilter(OutputFilter.SAM_SINGLE_END).matedMaxMismatches(IntegerOrPercentage.valueOf("100%")).create()).create();
    final NgsParams params = NgsParams.builder()
    .buildFirstParams(SequenceParams.builder().directory(inDir).useMemReader(true).create())
    .searchParams(SequenceParams.builder().directory(template).useMemReader(true).loadNames(true).create())
    .maxFragmentLength(1000).minFragmentLength(0)
    .outputParams(op)
    .gapOpenPenalty(1).gapExtendPenalty(1).substitutionPenalty(1).unknownsPenalty(0).alignerBandWidthFactor(new MaxShiftFactor(1))
    .create();
    final ByteArrayOutputStream log = new ByteArrayOutputStream();
    final PrintStream prLog = new PrintStream(log);
    Diagnostic.setLogStream(prLog);


    final MapQScoringReadBlocker blocker = new MapQScoringReadBlocker((int) params.buildFirstParams().reader().numberSequences(), 10);

    try (ByteArrayOutputStream outStream = new ByteArrayOutputStream()) {
      final MapStatistics stats = new SingleEndMapStatistics(null);
      final MyStatusReadIdListenerImpl listener = new MyStatusReadIdListenerImpl(2099, stats);
      try (SingleEndTempFileWriter w = new SingleEndTempFileWriter(params, listener, SharedResources.generateSharedResources(params))) {
        final HashingRegion region = new HashingRegion(0, 0, 0, 100, -1, -1);
        w.setClipRegion(region);

        w.initialiseAlignments(outStream, blocker);
//        w.initialiseUnmapped(outUnmappedStream, false);
        w.nextTemplateId(0);
        for (int i = 0; i < 2050; ++i) {
          w.alignmentResultUnfiltered(i, false, 1);
        }
        for (int i = 2050; i < 2099; ++i) {
          w.alignmentResultUnfiltered(i, false, 0);
        }
        listener.calculateStatistics(false, false);
//        TestUtils.containsAll(StringUtils.toString(stats), "UNMAPPED_READS=0", "MISSING_LEFT=0", "UNMAPPED_READS_PERCENT=0",
//            "MATED_UNIQUE_READS=0", "MISSING_RIGHT=0", "UNMATED_AMBIG_READS_PERCENT=100", "UNMATED_AMBIG_READS=2099", "MATED_UNIQUE_READS_PERCENT=0",
//            "TOTAL_READS=2099", "MAPPED_READS_PERCENT=100", "MAPPED_READS=2099");
        assertEquals(0L, stats.value(MapStatisticsField.UNMAPPED_NO_HITS, Arm.LEFT));
        assertEquals(0L, stats.value(MapStatisticsField.MISSING, Arm.LEFT));
        assertEquals(0L, stats.value(MapStatisticsField.MATED_UNIQUE_READS, Arm.LEFT));
        assertEquals(0L, stats.value(MapStatisticsField.MATED_AMBIG_READS, Arm.LEFT));
        assertEquals(2099L, stats.value(MapStatisticsField.UNMATED_AMBIG_READS, Arm.LEFT));
        assertEquals(2099L, stats.value(MapStatisticsField.TOTAL_READS, Arm.LEFT));
        assertEquals(0.0, stats.valueAsPercent(MapStatisticsField.UNMATED_UNIQUE_READS, Arm.LEFT));
        assertEquals(0.0, stats.valueAsPercent(MapStatisticsField.MISSING, Arm.LEFT));
        assertEquals(0.0, stats.valueAsPercent(MapStatisticsField.MATED_UNIQUE_READS, Arm.LEFT));
        assertEquals(0.0, stats.valueAsPercent(MapStatisticsField.MATED_AMBIG_READS, Arm.LEFT));
        assertEquals(100.0, stats.valueAsPercent(MapStatisticsField.UNMATED_AMBIG_READS, Arm.LEFT));
        assertEquals(100.0, stats.valueAsPercent(MapStatisticsField.TOTAL_READS, Arm.LEFT));
        assertNotNull(w.getBlocker());
      }
    } finally {
      Diagnostic.setLogStream();
      prLog.close();

    }

    TestUtils.containsAll(log.toString(),
        "Reordering buffer used capacity of 2099 records",
        "Statistics of multithreaded blocked",
        "2099 reads had count 0",
        "Total reads 2099",
        "Alignment score filter ",
        ": 2099 passed, 2099 total",
        "Duplicates detected during SAM writing: 0",
        "Freeing: 66 bytes",
        "AbstractSamAlignmentWriter.close() "
   //     "AbstractSamAlignmentWriter Allocating (whole): 66 bytes, startpos= 0 endpos= 66 newlength=66"
        );
  }
}
