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

import com.rtg.AbstractTest;
import com.rtg.index.hash.ngs.OutputProcessor;
import com.rtg.launcher.HashingRegion;
import com.rtg.launcher.ISequenceParams;
import com.rtg.launcher.MockReaderParams;
import com.rtg.launcher.MockSequenceParams;
import com.rtg.launcher.SequenceParams;
import com.rtg.mode.ProteinScoringMatrix;
import com.rtg.mode.SequenceMode;
import com.rtg.protein.ProteinOutputProcessor;
import com.rtg.protein.TopEqualProteinOutputProcessor;
import com.rtg.protein.TopEqualProteinOutputProcessorTest;
import com.rtg.protein.TopNProteinOutputProcessor;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.reader.SequencesReader;
import com.rtg.util.InvalidParamsException;
import com.rtg.util.TestUtils;
import com.rtg.util.io.TestDirectory;

/**
 *
 */
public class OutputFilterTest extends AbstractTest {

  private static final String DNA_FASTA = ">test\nacgtacgt";
  private static final String PROTEIN_FASTA = ">x\nacgt\n";

  public static ISequenceParams getSequenceParamsDna(String fasta) throws IOException {
    final SequencesReader dr = ReaderTestUtils.getReaderDnaMemory(fasta);
    return new MockSequenceParams(new MockReaderParams(dr), SequenceMode.BIDIRECTIONAL, 0, dr.numberSequences());
  }

  public static ISequenceParams getSequenceParamsProtein(String fasta) throws IOException {
    final SequencesReader pr = ReaderTestUtils.getReaderProteinMemory(fasta);
    return new MockSequenceParams(new MockReaderParams(pr), SequenceMode.PROTEIN, 0, pr.numberSequences());
  }


  public void test() {
    TestUtils.testPseudoEnum(OutputFilter.class, "[NONE, PAIRED_END, TOPN_PAIRED_END, PROTEIN_ALL_HITS, PROTEIN_TOPEQUAL, PROTEIN_TOPN, SAM_SINGLE_END, SAM_UNFILTERED, NULL]");
  }

  public void testNone() {
    final OutputFilter f = OutputFilter.NONE;
    assertEquals(OutputFilter.NONE, f);
    assertEquals(0, f.ordinal());
    assertEquals("NONE", f.toString());
    assertEquals(OutputFilter.NONE, OutputFilter.valueOf("NONE"));
  }

  public void testNull() {
    final OutputFilter f = OutputFilter.NULL;
    assertEquals(OutputFilter.NULL, f);
    assertEquals(8, f.ordinal());
    assertEquals("NULL", f.toString());
    assertEquals(OutputFilter.NULL, OutputFilter.valueOf("NULL"));
  }

  public void testNull1() throws IOException {
    final ISequenceParams sp = getSequenceParamsDna(PROTEIN_FASTA);
    final NgsParams ngsp = NgsParams.builder().buildFirstParams(sp).searchParams(sp).outputParams(new NgsTestUtils.OverriddenNgsOutputParams(NgsTestUtils.OverriddenNgsOutputParams.builder())).create();
    final OutputProcessor filter = OutputFilter.NULL.makeProcessor(ngsp, null);
    assertNotNull(filter);
    assertTrue(filter instanceof NullOutputProcessor);
  }

  public void testNone1() throws IOException {
    final ISequenceParams sp = getSequenceParamsDna(PROTEIN_FASTA);
    final NgsParams ngsp = NgsParams.builder().buildFirstParams(sp).searchParams(sp).outputParams(new NgsTestUtils.OverriddenNgsOutputParams(NgsTestUtils.OverriddenNgsOutputParams.builder())).create();
    final OutputProcessor filter = OutputFilter.NONE.makeProcessor(ngsp, null);
    assertNotNull(filter);
    assertTrue(filter instanceof DefaultOutputProcessorSynch);
  }

  public void testNone1Syn() throws IOException {
    final ISequenceParams sp = getSequenceParamsDna(PROTEIN_FASTA);
    final NgsParams ngsp = NgsParams.builder().numberThreads(2).buildFirstParams(sp).searchParams(sp).outputParams(new NgsTestUtils.OverriddenNgsOutputParams(NgsTestUtils.OverriddenNgsOutputParams.builder())).create();
    final OutputProcessor filter = OutputFilter.NONE.makeProcessor(ngsp, null);
    filter.close();
    assertNotNull(filter);
    assertTrue(filter instanceof DefaultOutputProcessorSynch);
  }


  public void testPairedEnd() throws IOException {
    try (final TestDirectory tdir = new TestDirectory("outputfilter")) {
      final File out = new File(tdir, "out");
      final File seq = ReaderTestUtils.getDNADir(">0\naacatcatcatcatactaacgt", tdir);
      final NgsFilterParams filterParams = NgsFilterParams.builder().outputFilter(OutputFilter.PAIRED_END).create();
      final NgsOutputParams ngsop = NgsOutputParams.builder().progress(false).outputDir(out).filterParams(filterParams).create();
      final SequenceParams param = SequenceParams.builder().directory(seq).useMemReader(true).loadNames(true).create();
      try (NgsParams ngsp = NgsParams.builder().numberThreads(1).buildFirstParams(param).searchParams(param).buildSecondParams(param).outputParams(ngsop).maskParams("amask").maxFragmentLength(1000).minFragmentLength(10).create()) {
        try (TopNPairedEndOutputProcessorSync topeq = (TopNPairedEndOutputProcessorSync) ngsop.outFilter().makeProcessor(ngsp, null)) {
          topeq.threadClone(HashingRegion.NONE).threadFinish();
          assertNotNull(topeq);
          topeq.finish();
        }
        assertFalse(new File(out, "unmated.sam").exists());
        assertFalse(new File(out, "unmapped.sam").exists());

        try (NgsParams ngsp2 = NgsParams.builder().numberThreads(2).buildFirstParams(param).searchParams(param).buildSecondParams(param).outputParams(ngsop).maxFragmentLength(1000).minFragmentLength(10).create()) {
          try (OutputProcessor op = ngsop.outFilter().makeProcessor(ngsp2, null)) {
            op.threadClone(HashingRegion.NONE).threadFinish();
            assertTrue(op instanceof TopNPairedEndOutputProcessorSync);
            op.finish();
          }
          assertFalse(new File(out, "unmated.sam").exists());
          assertFalse(new File(out, "unmapped.sam").exists());
        }
      } catch (final RuntimeException e) {
        fail(e.getMessage());
      }
    }

  }

  public void testTopNPairedEnd() throws Exception {
    try (final TestDirectory tdir = new TestDirectory("outputfilter")) {
      final File out = new File(tdir, "out");
      final File seq = ReaderTestUtils.getDNADir(">0\naacatcatcatcatactaacgt", tdir);
      final NgsFilterParams filterParams = NgsFilterParams.builder().outputFilter(OutputFilter.TOPN_PAIRED_END).create();
      final NgsOutputParams ngsop = NgsOutputParams.builder().progress(false).outputDir(out).filterParams(filterParams).create();
      final SequenceParams param = SequenceParams.builder().directory(seq).useMemReader(true).loadNames(true).create();
      try (NgsParams ngsp = NgsParams.builder().numberThreads(1).buildFirstParams(param).searchParams(param).buildSecondParams(param).outputParams(ngsop).maskParams("amask").maxFragmentLength(1000).minFragmentLength(10).create()) {
        final TopNPairedEndOutputProcessorSync topnPE = (TopNPairedEndOutputProcessorSync) ngsop.outFilter().makeProcessor(ngsp, null);
        topnPE.threadClone(HashingRegion.NONE).threadFinish();
        try {
          assertNotNull(topnPE);
          topnPE.finish();
        } finally {
          topnPE.close();
        }
      }
      try (NgsParams ngsp2 = NgsParams.builder().numberThreads(2).buildFirstParams(param).searchParams(param).buildSecondParams(param).outputParams(ngsop).maxFragmentLength(1000).minFragmentLength(10).create()) {
        final OutputProcessor op = ngsop.outFilter().makeProcessor(ngsp2, null);
        op.threadClone(HashingRegion.NONE).threadFinish();
        assertTrue(op instanceof TopNPairedEndOutputProcessorSync);
        op.finish();
        op.close();
      }
    }
  }

  public void testProteinOP() throws IOException, InvalidParamsException {
    try (final TestDirectory tmp = new TestDirectory("proteinop")) {
      final NgsParams params = NgsParams.builder().outputParams(NgsOutputParams.builder().outputDir(tmp).create())
          .buildFirstParams(getSequenceParamsDna(DNA_FASTA))
          .searchParams(getSequenceParamsProtein(PROTEIN_FASTA))
          .proteinScoringMatrix(new ProteinScoringMatrix())
          .create();
      final OutputProcessor p = OutputFilter.PROTEIN_ALL_HITS.makeProcessor(params, null);
      assertTrue(p instanceof ProteinOutputProcessor);
      p.finish();
      p.close();
      assertTrue(new File(tmp, "alignments.tsv").exists());
    }
  }

  public void testProteinTopEqual1() throws IOException, InvalidParamsException {
    try (final TestDirectory tmp = new TestDirectory("proteinop")) {
      final NgsParams params = NgsParams.builder().outputParams(NgsOutputParams.builder().outputDir(tmp).create())
          .buildFirstParams(getSequenceParamsDna(DNA_FASTA))
          .searchParams(getSequenceParamsProtein(PROTEIN_FASTA))
          .proteinScoringMatrix(new ProteinScoringMatrix())
          .create();
      final OutputProcessor p = OutputFilter.PROTEIN_TOPEQUAL.makeProcessor(params, null);
      assertTrue(p instanceof TopEqualProteinOutputProcessor);
      p.finish();
      p.close();
      assertTrue(new File(tmp, "alignments.tsv").exists());
    }
  }

  public void testProteinTopN1() throws IOException, InvalidParamsException {
    try (final TestDirectory tmp = new TestDirectory("proteinop")) {
      final NgsParams params = TopEqualProteinOutputProcessorTest.createParams(tmp, 5);
      final OutputProcessor p = OutputFilter.PROTEIN_TOPN.makeProcessor(params, null);
      assertTrue(p instanceof TopNProteinOutputProcessor);
      p.finish();
      p.close();
      assertTrue(new File(tmp, "alignments.tsv").exists());
    }

    try (final TestDirectory tmp2 = new TestDirectory("proteintopn")) {
      final NgsParams params = TopEqualProteinOutputProcessorTest.createParams(tmp2, 1);
      final OutputProcessor p = OutputFilter.PROTEIN_TOPN.makeProcessor(params, null);
      assertTrue(p instanceof TopEqualProteinOutputProcessor);
      p.finish();
      p.close();
      assertTrue(new File(tmp2, "alignments.tsv").exists());
    }
  }


  public void testProtein1() throws IOException, InvalidParamsException {
    try (TestDirectory outFileDir = new TestDirectory("outputFilter")) {
      final NgsFilterParams filterParams = NgsFilterParams.builder().outputFilter(OutputFilter.PROTEIN_ALL_HITS).create();
      final NgsOutputParams ngsop = NgsOutputParams.builder().progress(false).outputDir(outFileDir).filterParams(filterParams).create();
      final NgsParams ngsp = NgsParams.builder()
        .buildFirstParams(getSequenceParamsDna(DNA_FASTA))
        .searchParams(getSequenceParamsProtein(PROTEIN_FASTA))
        .outputParams(ngsop).proteinScoringMatrix(new ProteinScoringMatrix()).numberThreads(1).create();
      final OutputProcessor prot = ngsop.outFilter().makeProcessor(ngsp, null);
      prot.close();
      assertNotNull(prot);
      assertTrue(prot instanceof ProteinOutputProcessor);
    }
  }

  public void testSingle1() throws IOException, InvalidParamsException {
    try (final TestDirectory tmp = new TestDirectory("unfilteredope")) {
      final NgsParams params = NgsParams.builder().outputParams(NgsOutputParams.builder().outputDir(tmp).create())
          .buildFirstParams(getSequenceParamsDna(DNA_FASTA))
          .searchParams(getSequenceParamsDna(DNA_FASTA))
          .create();
      final OutputProcessor p = OutputFilter.SAM_SINGLE_END.makeProcessor(params, null);
      p.threadClone(HashingRegion.NONE).threadFinish();
      assertTrue(p instanceof SamSingleEndOutputProcessor);
      p.finish();
      p.close();
      assertTrue(new File(tmp, "alignments.sam").exists());
    }
  }

  public void testUnfiltered1() throws IOException, InvalidParamsException {
    try (final TestDirectory tmp = new TestDirectory("unfilteredope")) {
      final NgsParams params = NgsParams.builder().outputParams(NgsOutputParams.builder().outputDir(tmp).create())
          .buildFirstParams(getSequenceParamsDna(DNA_FASTA))
          .searchParams(getSequenceParamsDna(DNA_FASTA))
          .create();
      final OutputProcessor p = OutputFilter.SAM_UNFILTERED.makeProcessor(params, null);
      p.threadClone(HashingRegion.NONE).threadFinish();
      assertTrue(p instanceof UnfilteredSingleEndOutputProcessor);
      p.finish();
      p.close();
      assertTrue(new File(tmp, "alignments.sam").exists());
    }
  }

  public void testUnfiltered2() throws IOException, InvalidParamsException {
    try (final TestDirectory tmp = new TestDirectory("unfilteredse")) {
      try (NgsParams params = NgsParams.builder().outputParams(NgsOutputParams.builder().outputDir(tmp).create())
        .buildFirstParams(getSequenceParamsDna(DNA_FASTA))
        .buildSecondParams(getSequenceParamsDna(DNA_FASTA))
        .searchParams(getSequenceParamsDna(DNA_FASTA))
        .maxFragmentLength(500).minFragmentLength(100)
        .create()) {
        try (OutputProcessor p = OutputFilter.SAM_UNFILTERED.makeProcessor(params, null)) {
          p.threadClone(HashingRegion.NONE).threadFinish();
          assertTrue(p instanceof UnfilteredPairedEndOutputProcessor);
          p.finish();
        }
      }
      assertTrue(new File(tmp, "alignments.sam").exists());
    }
  }

  public void testHashEquals() {
    TestUtils.equalsHashTest(new OutputFilter[][] {{OutputFilter.NONE}, {OutputFilter.PAIRED_END}, {OutputFilter.TOPN_PAIRED_END}});
  }

}
