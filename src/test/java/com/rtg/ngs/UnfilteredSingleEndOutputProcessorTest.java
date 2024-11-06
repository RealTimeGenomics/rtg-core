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


import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.PrintStream;

import com.rtg.AbstractTest;
import com.rtg.index.hash.ngs.OutputProcessor;
import com.rtg.launcher.HashingRegion;
import com.rtg.launcher.SequenceParams;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.util.StringUtils;
import com.rtg.util.TestUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.IOUtils;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.util.test.FileHelper;

/**
 * Tests corresponding class
 */
public class UnfilteredSingleEndOutputProcessorTest extends AbstractTest {

  protected File mDir = null;
  @Override
  public void setUp() throws IOException {
    super.setUp();
    mDir = FileHelper.createTempDirectory();
  }

  @Override
  public void tearDown() throws IOException {
    super.tearDown();
    assertTrue(!mDir.exists() || FileHelper.deleteAll(mDir));
    mDir = null;
  }

  static final String TEMP = ">t" + StringUtils.LS + "tgcaagacaagagggcctcc" + StringUtils.LS;
  static final String READ1 = "tgcaagacaagagggcctcc";
  static final String READ2 = "ggaggccctcttgtcttgca";
  static final String READS = ">r1" + StringUtils.LS + READ1 + StringUtils.LS
  + ">r2" + StringUtils.LS + READ2 + StringUtils.LS;

  protected NgsParamsBuilder getDefaultBuilder(String reads) throws IOException {
    return getDefaultBuilder(null, false, reads);
  }

  private NgsParamsBuilder getDefaultBuilder(final File tempDir, final boolean gzipOutputs, String reads) throws IOException {
    return getDefaultBuilder(false, tempDir, gzipOutputs, reads);
  }

  private NgsParamsBuilder getDefaultBuilder(final boolean sdf, final File tempDir, final boolean gzipOutputs, String reads) throws IOException {
    final File templateok = FileUtils.createTempDir("template", "ngs", mDir);
    final File readsok = FileUtils.createTempDir("reads", "ngs", mDir);
    final File hitsDir = new File(mDir, "hitDir");

    ReaderTestUtils.getReaderDNA(TEMP, templateok, null).close();
    ReaderTestUtils.getReaderDNA(reads, readsok, null).close();

    final NgsFilterParams filterParams = NgsFilterParams.builder().outputFilter(OutputFilter.PAIRED_END)
    .topN(10).errorLimit(5).zip(gzipOutputs).create();
    final NgsOutputParams outputParams = NgsOutputParams.builder()
    .progress(false).outputDir(hitsDir).sam(!sdf).sdf(sdf).outputReadNames(true)
    .tempFilesDir(tempDir).filterParams(filterParams).create();

    return NgsParams.builder()
    .buildFirstParams(SequenceParams.builder().directory(readsok).useMemReader(true).loadNames(true).loadFullNames(true).create())
    .searchParams(SequenceParams.builder().directory(templateok).useMemReader(true).loadNames(true).create())
    .outputParams(outputParams)
    .substitutionPenalty(1).gapOpenPenalty(1).gapExtendPenalty(1).unknownsPenalty(0)
    .maxFragmentLength(1000).minFragmentLength(0);
  }

  //  static class DummyPairedEndOutputProcessor extends PairedEndOutputProcessor {
  //
  //    public DummyPairedEndOutputProcessor(final NgsParams param) {
  //      super(param, null, param.outStream());
  //    }
  //  }

  public void testStandard() throws Exception {
    final ByteArrayOutputStream log = new ByteArrayOutputStream();
    try (NgsParams param = getDefaultBuilder(READS).create()) {
      final PrintStream prLog = new PrintStream(log);
      Diagnostic.setLogStream(prLog);
      try {
        final UnfilteredSingleEndOutputProcessor suop = new UnfilteredSingleEndOutputProcessor(param, null, true);
        try {
          TestUtils.containsAll(suop.toString(), "SamUnfilteredOutputProcessor");
          suop.process(0, "F", 0, 1, 0, 0);
          //sseop.process(0, "R", 1, 16, 0, 0);

          suop.finish();
          suop.close();
          ///assertNull(peop.mSamWriter);
          final String outStr = IOUtils.readAll(new File(param.outputParams().directory(), "alignments.sam"));
          assertTrue(outStr, outStr.contains("r1\t0\tt\t1\t255\t20=\t*\t0\t0\tTGCAAGACAAGAGGGCCTCC\t*\tAS:i:0\tNM:i:0"));
          assertTrue(new File(param.outputParams().directory(), "unmapped.sam").isFile());
          suop.mUnmappedTracker.calculateStatistics(param.paired(), false);
          prLog.flush();

          //        TestUtils.containsAll(baos.toString(), new String[] {
          //          "Sliding window collector statistics",
          //          "Setting stats:",
          //          "hits = 2"
          //        });
        } finally {
          suop.close();
        }
      } finally {
        prLog.close();
        //System.err.println("log:\n" + log.toString());
      }
    }
//    TestUtils.containsAll(log.toString(), new String[]{"Extracting hits for reads",
//        "AlignmentOutput",
//        "Processing ",
//        " hits for reads",
//      "Merging alignment results"});
//      Diagnostic.setLogStream();
  }

  public void testNoSamSDF() throws Exception {
    final ByteArrayOutputStream log = new ByteArrayOutputStream();
    try (NgsParams param = getDefaultBuilder(true, null, false, READS).create()) {
      final PrintStream prLog = new PrintStream(log);
      Diagnostic.setLogStream(prLog);
      try {
        final UnfilteredSingleEndOutputProcessor suop = new UnfilteredSingleEndOutputProcessor(param, null, true);
        try {
          TestUtils.containsAll(suop.toString(), "SamUnfilteredOutputProcessor");
          suop.process(0, "F", 0, 1, 0, 0);

          suop.finish();
          suop.close();
          assertFalse(new File(param.outputParams().directory(), "alignments.sam").isFile());
          assertFalse(new File(param.outputParams().directory(), "unmapped.sam").isFile());
          assertTrue(new File(param.outputParams().directory(), "alignments.sdf").isDirectory());
          assertTrue(new File(param.outputParams().directory(), "unmapped.sdf").isDirectory());
          assertEquals(0, param.outputParams().directory().list(new FilenameFilter() {
            @Override
            public boolean accept(File dir, String name) {
              return name.contains("TEMP_SAM");
            }
          }).length);
          prLog.flush();

        } finally {
          suop.close();
        }
      } finally {
        prLog.close();
      }
    }
  }

  public void testInnerUnfilteredOutputProcessor() throws Exception {
    try (MemoryPrintStream log = new MemoryPrintStream()) {
      Diagnostic.setLogStream(log.printStream());
      final NgsParams param = getDefaultBuilder(READS).create();
      final UnfilteredSingleEndOutputProcessor suop = new UnfilteredSingleEndOutputProcessor(param, null, false);
      final OutputProcessor inner = suop.threadClone(HashingRegion.NONE);

      try {
        inner.threadClone(HashingRegion.NONE);
        fail();
      } catch (final UnsupportedOperationException e) {
        assertEquals("Should never get called.", e.getMessage());
      }
      inner.threadFinish();
      assertTrue(log.toString().contains("Child finish"));
    } finally {
      Diagnostic.setLogStream();

    }
  }

}
