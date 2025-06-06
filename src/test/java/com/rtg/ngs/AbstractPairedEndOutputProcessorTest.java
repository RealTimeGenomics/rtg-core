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

import com.rtg.launcher.AbstractNanoTest;
import com.rtg.launcher.SequenceParams;
import com.rtg.ngs.tempstage.BinaryTempFileRecord;
import com.rtg.ngs.tempstage.PairedTempFileWriterImpl;
import com.rtg.ngs.tempstage.TempRecordReader;
import com.rtg.ngs.tempstage.TempRecordReaderNio;
import com.rtg.pairedend.SlidingWindowCollector;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.sam.SamCommandHelper;
import com.rtg.util.InvalidParamsException;
import com.rtg.util.StringUtils;
import com.rtg.util.TestUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.util.machine.MachineOrientation;
import com.rtg.util.test.FileHelper;

/**
 */
public abstract class AbstractPairedEndOutputProcessorTest extends AbstractNanoTest {

  protected File mDir = null;

  @Override
  public void setUp() throws IOException {
    super.setUp();
    mDir = FileHelper.createTempDirectory();
  }

  @Override
  public void tearDown() throws IOException {
    assertTrue(mDir == null || !mDir.exists() || FileHelper.deleteAll(mDir));
    mDir = null;
    super.tearDown();
  }

  static final String TEMPLATE = ">t" + StringUtils.LS + "tgcaagacaagagggcctcc" + StringUtils.LS;
  static final String TEMP_LEFT = "tgcaagacaagagggcctcc";
  static final String TEMP_RIGHT = "ggaggccctcttgtcttgca";
  static final String READ_LEFT = ">r" + StringUtils.LS + TEMP_LEFT + StringUtils.LS;
  static final String READ_RIGHT = ">r" + StringUtils.LS + TEMP_RIGHT + StringUtils.LS;


  protected NgsParamsBuilder getDefaultBuilder() throws IOException, InvalidParamsException {
    return getDefaultBuilder(null, false, OutputFilter.TOPN_PAIRED_END, null);
  }

  NgsParamsBuilder getDefaultBuilder(File tempDir, boolean gzipOutputs, OutputFilter outputFilter, File readGroupFile) throws IOException, InvalidParamsException {
    final File templateok = FileUtils.createTempDir("template", "ngs", mDir);
    final File leftok = FileUtils.createTempDir("left", "ngs", mDir);
    final File rightok = FileUtils.createTempDir("right", "ngs", mDir);
    final File hitsDir = new File(mDir, "hitDir");

    ReaderTestUtils.getReaderDNA(TEMPLATE, templateok, null).close();
    ReaderTestUtils.getReaderDNA(READ_LEFT, leftok, null).close();
    ReaderTestUtils.getReaderDNA(READ_RIGHT, rightok, null).close();

    final NgsFilterParams filterParams = NgsFilterParams.builder().outputFilter(outputFilter)
        .topN(10).errorLimit(5).zip(gzipOutputs).create();
    final NgsOutputParamsBuilder outputParamsBuilder = NgsOutputParams.builder()
        .progress(false).outputDir(hitsDir)
        .tempFilesDir(tempDir).filterParams(filterParams);
    if (readGroupFile == null) {
      outputParamsBuilder.calibrate(false);
    } else {
      outputParamsBuilder.calibrate(true);
      outputParamsBuilder.readGroup(SamCommandHelper.validateAndCreateSamRG(readGroupFile.toString(), SamCommandHelper.ReadGroupStrictness.REQUIRED));
    }
    final NgsOutputParams outputParams = outputParamsBuilder.create();
    return NgsParams.builder()
        .buildFirstParams(SequenceParams.builder().directory(leftok).useMemReader(true).create())
        .buildSecondParams(SequenceParams.builder().directory(rightok).useMemReader(true).create())
        .searchParams(SequenceParams.builder().directory(templateok).useMemReader(true).loadNames(true).create())
        .outputParams(outputParams)
        .maskParams(new NgsMaskParamsGeneral(4, 1, 1, 1))
        .substitutionPenalty(1).gapOpenPenalty(1).gapExtendPenalty(1).unknownsPenalty(0)
        .maxFragmentLength(1000).minFragmentLength(0);
  }

  PairedEndOutputProcessor getOutputProc(NgsParams param, File tempFile) throws IOException {
    final ReadStatusTracker mUnmappedTracker = new ReadStatusTracker((int) param.buildFirstParams().numberSequences(), null);

    final SharedResources sharedResources = SharedResources.generateSharedResources(param);

    final PairedTempFileWriterImpl mSamWriter = new PairedTempFileWriterImpl(param, mUnmappedTracker, sharedResources);
    assertTrue(param.outputParams().directory().mkdir());
    mSamWriter.initialiseMated(FileUtils.createOutputStream(tempFile, param.compressOutput()));
    final int max = param.maxFragmentLength();
    final int min = param.minFragmentLength();

    final SlidingWindowCollector collector = new SlidingWindowCollector(max, min, MachineOrientation.ANY, mSamWriter, sharedResources, param.outputParams().calibrateRegions());

    return getPEOP(mSamWriter, collector);
  }

  abstract PairedEndOutputProcessor getPEOP(PairedTempFileWriterImpl writer, SlidingWindowCollector collector);

  public void testStuff() throws Exception {
    try (NgsParams param = getDefaultBuilder().create()) {
      final MemoryPrintStream mps = new MemoryPrintStream();
      Diagnostic.setLogStream(mps.printStream());
      final File tempFile = new File(param.outputParams().directory(), "mated.sam");
      final PairedEndOutputProcessor peop = getOutputProc(param, tempFile);
      try {
        assertNotNull(peop.mSamWriter);

        peop.process(0, false, 0, 1);
        peop.process(0, true, 1, 3);

        peop.finish();
      } finally {
        peop.close();
      }
      assertNull(peop.mSamWriter);

      try (TempRecordReader trr = new TempRecordReaderNio(FileUtils.createInputStream(tempFile, false), new TempRecordReader.RecordFactory(true, false, false, false))) {
        BinaryTempFileRecord bar = trr.readRecord();
        assertNotNull(bar);
        assertEquals(0, bar.getReadId());
        assertEquals(147, bar.getSamFlags() & 0xff);
        assertEquals(0, bar.getReferenceId());
        assertEquals(1, bar.getStartPosition());
        assertEquals(-20, bar.getTemplateLength());

        bar = trr.readRecord();
        assertNotNull(bar);
        assertEquals(0, bar.getReadId());
        assertEquals(99, bar.getSamFlags() & 0xff);
        assertEquals(0, bar.getReferenceId());
        assertEquals(1, bar.getStartPosition());
        assertEquals(20, bar.getTemplateLength());

        assertNull(trr.readRecord());
      }


      //0 147 0 1 255 20= = 1 -20 * * AS:i:0  NM:i:0  MQ:i:255  XA:i:0
      //0 99  0 1 255 20= = 1 20  * * AS:i:0  NM:i:0  MQ:i:255  XA:i:0


      TestUtils.containsAll(mps.toString(), "Sliding window collector statistics",
        "hits = 2");
    }
  }

}
