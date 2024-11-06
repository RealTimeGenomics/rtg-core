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

import com.rtg.index.hash.ngs.OutputProcessor;
import com.rtg.launcher.HashingRegion;
import com.rtg.launcher.SequenceParams;
import com.rtg.launcher.globals.GlobalFlags;
import com.rtg.ngs.blocking.MapQScoringReadBlocker;
import com.rtg.reader.IndexFile;
import com.rtg.reader.NamesInterface;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.reader.Sdf2Fasta;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.util.test.FileHelper;

import junit.framework.TestCase;

/**
 */
public class DummySdfOutputProcessorTest extends TestCase {

  private class BlahOutputProcessor extends AbstractSdfOutputProcessor {

    BlahOutputProcessor(NgsParams param) throws IOException {
      super(param, null, false, false);
    }

    @Override
    public void process(long templateName, String frame, int readId, int tStart, int score, int scoreIndel) {

    }

    @Override
    public OutputProcessor threadClone(HashingRegion region) {
      return null;
    }

    @Override
    public void threadFinish() {
    }

    @Override
    public void close() {
    }

    @Override
    protected FilterConcatIntermediateFiles filterConcatNonMated(MapQScoringReadBlocker blockerLeft, MapQScoringReadBlocker blockerRight, File[] tempFiles, SingleEndTopRandomImplementation.HitRecord[] hitsToKeep, NamesInterface templateNames, File outFile) {
      throw new UnsupportedOperationException();
    }
  }

  @Override
  protected void tearDown() throws Exception {
    super.tearDown();
    GlobalFlags.resetAccessedStatus();
  }

  static final String TEMPLATE = ">t\n" + "tgcaagacaagagggcctcc\n";
  static final String TEMP_LEFT = "TGCAAGACAAGAGGGCCTCC";
  static final String TEMP_RIGHT = "GGAGGCCCTCTTGTCTTGCA";
  static final String READ_LEFT = ">r\n" + TEMP_LEFT + "\n"
                                + ">r1\n" + TEMP_LEFT + "\n"
                                + ">r2\n" + TEMP_LEFT + "\n";
  static final String READ_RIGHT = ">r\n" + TEMP_RIGHT + "\n"
                                 + ">r1\n" + TEMP_RIGHT + "\n"
                                 + ">r2\n" + TEMP_RIGHT + "\n";

  public void testSDFOutput() throws Exception {
    final File tmpDir = FileUtils.createTempDir("tmp", "dir");
    try {
      final File templateok = FileUtils.createTempDir("template", "ngs", tmpDir);
      final File leftok = FileUtils.createTempDir("left", "ngs", tmpDir);
      final File rightok = FileUtils.createTempDir("right", "ngs", tmpDir);
      final File hitsDir = new File(tmpDir, "hitDir");

      ReaderTestUtils.getReaderDNA(TEMPLATE, templateok, null).close();
      ReaderTestUtils.getReaderDNA(READ_LEFT, leftok, null).close();
      ReaderTestUtils.getReaderDNA(READ_RIGHT, rightok, null).close();

      final NgsFilterParams filterParams = NgsFilterParams.builder().outputFilter(OutputFilter.NONE)
                                                                    .topN(10).errorLimit(5).zip(false).create();
      final NgsOutputParams outputParams = NgsOutputParams.builder()
        .progress(false).outputDir(hitsDir)
        .tempFilesDir(tmpDir).filterParams(filterParams).sdf(true).create();

      final NgsParams param = NgsParams.builder()
          .buildFirstParams(SequenceParams.builder().directory(leftok).useMemReader(true).loadNames(true).loadFullNames(true).create())
          .buildSecondParams(SequenceParams.builder().directory(rightok).useMemReader(true).loadNames(true).loadFullNames(true).create())
          .searchParams(SequenceParams.builder().directory(templateok).useMemReader(true).loadNames(true).create())
          .outputParams(outputParams)
          .maskParams(new NgsMaskParamsGeneral(4, 1, 1, 1))
        .maxFragmentLength(1000).minFragmentLength(0).create();

      final BlahOutputProcessor bop = new BlahOutputProcessor(param);

      bop.mUnmappedTracker.addStatus(0, ReadStatusTracker.UNMAPPED_FIRST);
      bop.mUnmappedTracker.addStatus(0, ReadStatusTracker.UNMAPPED_SECOND);
      bop.mUnmappedTracker.addStatus(1, ReadStatusTracker.UNIQUELY_MAPPED_FIRST);
      bop.mUnmappedTracker.addStatus(1, ReadStatusTracker.UNMAPPED_SECOND);

      bop.finish();

      final MemoryPrintStream mps = new MemoryPrintStream();
      final File alignmentsSdf = new File(hitsDir, "alignments.sdf");
      GlobalFlags.resetAccessedStatus();
      final int rc = new Sdf2Fasta().mainInit(new String[] {"-i", alignmentsSdf.getPath(), "-o", new File(hitsDir, "outfa").getPath(), "-Z"}, mps.outputStream(), mps.printStream());
      assertEquals(mps.toString(), 0, rc);

      assertEquals(">r1\n" + TEMP_LEFT + "\n"
                  + ">r2\n" + TEMP_LEFT + "\n", FileUtils.fileToString(new File(hitsDir, "outfa_1.fasta")));
      assertEquals(">r1\n" + TEMP_RIGHT + "\n"
                  + ">r2\n" + TEMP_RIGHT + "\n", FileUtils.fileToString(new File(hitsDir, "outfa_2.fasta")));

      assertEquals(mps.toString(), 0, new Sdf2Fasta().mainInit(new String[] {"-i", new File(hitsDir, "unmapped.sdf").getPath(), "-o", new File(hitsDir, "unmappedfa").getPath(), "-Z"}, mps.outputStream(), mps.printStream()));

      assertEquals(">r\n" + TEMP_LEFT + "\n", FileUtils.fileToString(new File(hitsDir, "unmappedfa_1.fasta")));
      assertEquals(">r\n" + TEMP_RIGHT + "\n", FileUtils.fileToString(new File(hitsDir, "unmappedfa_2.fasta")));

      //test the sdf ids are the same -_-
      mps.reset();

      final IndexFile ifl = new IndexFile(new File(alignmentsSdf, "left"));
      final IndexFile ifr = new IndexFile(new File(alignmentsSdf, "right"));

      assertEquals(ifl.getSdfId(), ifr.getSdfId());

    } finally {
      FileHelper.deleteAll(tmpDir);
    }
  }
}
