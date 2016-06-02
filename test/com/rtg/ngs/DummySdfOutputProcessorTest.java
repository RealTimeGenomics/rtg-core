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

import java.io.File;
import java.io.IOException;

import com.rtg.index.hash.ngs.OutputProcessor;
import com.rtg.launcher.globals.GlobalFlags;
import com.rtg.launcher.HashingRegion;
import com.rtg.launcher.SequenceParams;
import com.rtg.ngs.blocking.MapQScoringReadBlocker;
import com.rtg.reader.IndexFile;
import com.rtg.reader.PrereadNamesInterface;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.reader.Sdf2Fasta;
import com.rtg.util.StringUtils;
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
    protected FilterConcatIntermediateFiles filterConcatNonMated(MapQScoringReadBlocker blockerLeft, MapQScoringReadBlocker blockerRight, File[] tempFiles, SingleEndTopRandomImplementation.HitRecord[] hitsToKeep, PrereadNamesInterface templateNames, File outFile) {
      throw new UnsupportedOperationException();
    }
  }

  @Override
  protected void tearDown() throws Exception {
    super.tearDown();
    GlobalFlags.resetAccessedStatus();
  }

  static final String TEMPLATE = ">t" + StringUtils.LS + "tgcaagacaagagggcctcc" + StringUtils.LS;
  static final String TEMP_LEFT = "TGCAAGACAAGAGGGCCTCC";
  static final String TEMP_RIGHT = "GGAGGCCCTCTTGTCTTGCA";
  static final String READ_LEFT = ">r" + StringUtils.LS + TEMP_LEFT + StringUtils.LS
                                + ">r1" + StringUtils.LS + TEMP_LEFT + StringUtils.LS
                                + ">r2" + StringUtils.LS + TEMP_LEFT + StringUtils.LS;
  static final String READ_RIGHT = ">r" + StringUtils.LS + TEMP_RIGHT + StringUtils.LS
                                 + ">r1" + StringUtils.LS + TEMP_RIGHT + StringUtils.LS
                                 + ">r2" + StringUtils.LS + TEMP_RIGHT + StringUtils.LS;

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

      assertEquals(">r1" + StringUtils.LS + TEMP_LEFT + StringUtils.LS
                  + ">r2" + StringUtils.LS + TEMP_LEFT + StringUtils.LS, FileUtils.fileToString(new File(hitsDir, "outfa_1.fasta")));
      assertEquals(">r1" + StringUtils.LS + TEMP_RIGHT + StringUtils.LS
                  + ">r2" + StringUtils.LS + TEMP_RIGHT + StringUtils.LS, FileUtils.fileToString(new File(hitsDir, "outfa_2.fasta")));

      assertEquals(mps.toString(), 0, new Sdf2Fasta().mainInit(new String[] {"-i", new File(hitsDir, "unmapped.sdf").getPath(), "-o", new File(hitsDir, "unmappedfa").getPath(), "-Z"}, mps.outputStream(), mps.printStream()));

      assertEquals(">r" + StringUtils.LS + TEMP_LEFT + StringUtils.LS, FileUtils.fileToString(new File(hitsDir, "unmappedfa_1.fasta")));
      assertEquals(">r" + StringUtils.LS + TEMP_RIGHT + StringUtils.LS, FileUtils.fileToString(new File(hitsDir, "unmappedfa_2.fasta")));

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
