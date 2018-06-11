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

import static com.rtg.util.StringUtils.LS;

import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;

import com.rtg.index.hash.ngs.OutputProcessor;
import com.rtg.launcher.ISequenceParams;
import com.rtg.launcher.MockReaderParams;
import com.rtg.launcher.MockSequenceParams;
import com.rtg.launcher.ReaderParams;
import com.rtg.launcher.HashingRegion;
import com.rtg.launcher.SequenceParams;
import com.rtg.mode.DNAFastaSymbolTable;
import com.rtg.mode.ProteinScoringMatrix;
import com.rtg.mode.SequenceMode;
import com.rtg.protein.ProteinOutputProcessor;
import com.rtg.protein.ProteinOutputProcessorTest;
import com.rtg.protein.TopEqualProteinOutputProcessor;
import com.rtg.protein.TopEqualProteinOutputProcessorTest;
import com.rtg.protein.TopNProteinOutputProcessor;
import com.rtg.reader.CompressedMemorySequencesReader;
import com.rtg.reader.FastaSequenceDataSource;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.reader.SequencesReader;
import com.rtg.util.InvalidParamsException;
import com.rtg.util.TestUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.FileUtils;
import com.rtg.util.test.FileHelper;

import junit.framework.TestCase;

/**
 *
 */
public class OutputFilterTest extends TestCase {

  private static final byte[] FASTA = (">x" + LS + "acgt" + LS).getBytes();

  private SequencesReader getReaderDNA() throws IOException {
    final ArrayList<InputStream> inputStreams = new ArrayList<>();
    inputStreams.add(new ByteArrayInputStream(FASTA));
    final FastaSequenceDataSource ds = new FastaSequenceDataSource(inputStreams, new DNAFastaSymbolTable());
    return CompressedMemorySequencesReader.createSequencesReader(ds);
  }

  @Override
  public void setUp() {
    Diagnostic.setLogStream();
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
    final SequencesReader r = getReaderDNA();
    final ReaderParams rp = new MockReaderParams(r, SequenceMode.BIDIRECTIONAL);
    final ISequenceParams sp = new MockSequenceParams(rp, 0, r.numberSequences());
    final NgsParams ngsp = NgsParams.builder().buildFirstParams(sp).searchParams(sp).outputParams(new NgsTestUtils.OverriddenNgsOutputParams(NgsTestUtils.OverriddenNgsOutputParams.builder())).create();
    final OutputProcessor filter = OutputFilter.NULL.makeProcessor(ngsp, null);
    assertNotNull(filter);
    assertTrue(filter instanceof NullOutputProcessor);
  }

  public void testNone1() throws IOException {
    final SequencesReader r = getReaderDNA();
    final ReaderParams rp = new MockReaderParams(r, SequenceMode.BIDIRECTIONAL);
    final ISequenceParams sp = new MockSequenceParams(rp, 0, r.numberSequences());
    final NgsParams ngsp = NgsParams.builder().buildFirstParams(sp).searchParams(sp).outputParams(new NgsTestUtils.OverriddenNgsOutputParams(NgsTestUtils.OverriddenNgsOutputParams.builder())).create();
    final OutputProcessor filter = OutputFilter.NONE.makeProcessor(ngsp, null);
    assertNotNull(filter);
    assertTrue(filter instanceof DefaultOutputProcessorSynch);
  }

  public void testNone1Syn() throws IOException {
    final SequencesReader r = getReaderDNA();
    final ReaderParams rp = new MockReaderParams(r, SequenceMode.BIDIRECTIONAL);
    final ISequenceParams sp = new MockSequenceParams(rp, 0, r.numberSequences());
    final NgsParams ngsp = NgsParams.builder().numberThreads(2).buildFirstParams(sp).searchParams(sp).outputParams(new NgsTestUtils.OverriddenNgsOutputParams(NgsTestUtils.OverriddenNgsOutputParams.builder())).create();
    final OutputProcessor filter = OutputFilter.NONE.makeProcessor(ngsp, null);
    filter.close();
    assertNotNull(filter);
    assertTrue(filter instanceof DefaultOutputProcessorSynch);
  }


  public void testPairedEnd() throws IOException {
    final File seq = ReaderTestUtils.getDNADir(">0\naacatcatcatcatactaacgt");

    final NgsFilterParams filterParams = NgsFilterParams.builder().outputFilter(OutputFilter.PAIRED_END).create();
    final File out = FileUtils.createTempDir("pairedend", "outputfilter");
    try {
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
      } finally {
        assertTrue(FileHelper.deleteAll(seq));
      }
    } finally {
      assertTrue(FileHelper.deleteAll(out));
    }

  }

  public void testTopNPairedEnd() throws Exception {
    final File seq = ReaderTestUtils.getDNADir(">0\naacatcatcatcatactaacgt");

    final NgsFilterParams filterParams = NgsFilterParams.builder().outputFilter(OutputFilter.TOPN_PAIRED_END).create();
    final File out = FileUtils.createTempDir("pairedend", "outputfilter");
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

      try (NgsParams ngsp2 = NgsParams.builder().numberThreads(2).buildFirstParams(param).searchParams(param).buildSecondParams(param).outputParams(ngsop).maxFragmentLength(1000).minFragmentLength(10).create()) {
        final OutputProcessor op = ngsop.outFilter().makeProcessor(ngsp2, null);
        op.threadClone(HashingRegion.NONE).threadFinish();
        assertTrue(op instanceof TopNPairedEndOutputProcessorSync);
        op.finish();
        op.close();
      }
    } finally {
      assertTrue(FileHelper.deleteAll(seq));
      assertTrue(FileHelper.deleteAll(out));
    }

  }


  public void testProteinOP() throws IOException, InvalidParamsException {
    final File tmp = FileUtils.createTempDir("proteinop", "filter");
    try {
      final File input = new File(tmp, "1");
      ReaderTestUtils.getReaderDNA(">test\nacgtacgt", input, null).close();
      final File input2 = new File(tmp, "2");
      ReaderTestUtils.getReaderProtein(">test\nacgtacgt", input2).close();

      final NgsParams params = NgsParams.builder().outputParams(NgsOutputParams.builder().outputDir(tmp).create())
          .buildFirstParams(SequenceParams.builder().directory(input).useMemReader(true).create())
          .searchParams(SequenceParams.builder().directory(input2).loadNames(true).useMemReader(true).create())
          .proteinScoringMatrix(new ProteinScoringMatrix())
          .create();
      final OutputProcessor p = OutputFilter.PROTEIN_ALL_HITS.makeProcessor(params, null);
      assertTrue(p instanceof ProteinOutputProcessor);
      p.finish();
      p.close();
      assertTrue(new File(tmp, "alignments.tsv").exists());
    } finally {
      assertTrue(FileHelper.deleteAll(tmp));
    }
  }

  public void testProteinTopEqual1() throws IOException, InvalidParamsException {
    final File tmp = FileUtils.createTempDir("proteintope", "filter");
    try {
      final File input = new File(tmp, "1");
      ReaderTestUtils.getReaderDNA(">test\nacgtacgt", input, null).close();
      final File input2 = new File(tmp, "2");
      ReaderTestUtils.getReaderProtein(">test\nacgtacgt", input2).close();

      final NgsParams params = NgsParams.builder().outputParams(NgsOutputParams.builder().outputDir(tmp).create())
          .buildFirstParams(SequenceParams.builder().directory(input).useMemReader(true).create())
          .searchParams(SequenceParams.builder().loadNames(true).directory(input2).useMemReader(true).create())
          .proteinScoringMatrix(new ProteinScoringMatrix())
          .create();
      final OutputProcessor p = OutputFilter.PROTEIN_TOPEQUAL.makeProcessor(params, null);
      assertTrue(p instanceof TopEqualProteinOutputProcessor);
      p.finish();
      p.close();
      assertTrue(new File(tmp, "alignments.tsv").exists());
    } finally {
      assertTrue(FileHelper.deleteAll(tmp));
    }
  }

  public void testProteinTopN1() throws IOException, InvalidParamsException {
    final File tmp = FileUtils.createTempDir("proteintopn", "filter");
    try {
      final NgsParams params = TopEqualProteinOutputProcessorTest.createParams(tmp, ProteinOutputProcessorTest.READS_FASTA_PERFECT, ProteinOutputProcessorTest.TEMPLATE_FASTA, 5);
      final OutputProcessor p = OutputFilter.PROTEIN_TOPN.makeProcessor(params, null);
      assertTrue(p instanceof TopNProteinOutputProcessor);
      p.finish();
      p.close();
      assertTrue(new File(tmp, "alignments.tsv").exists());
    } finally {
      assertTrue(FileHelper.deleteAll(tmp));
    }

    final File tmp2 = FileUtils.createTempDir("proteintopn", "filter");
    try {
      final NgsParams params = TopEqualProteinOutputProcessorTest.createParams(tmp2, ProteinOutputProcessorTest.READS_FASTA_PERFECT, ProteinOutputProcessorTest.TEMPLATE_FASTA, 1);
      final OutputProcessor p = OutputFilter.PROTEIN_TOPN.makeProcessor(params, null);
      assertTrue(p instanceof TopEqualProteinOutputProcessor);
      p.finish();
      p.close();
      assertTrue(new File(tmp2, "alignments.tsv").exists());
    } finally {
      assertTrue(FileHelper.deleteAll(tmp2));
    }
  }


  public void testProtein1() throws IOException, InvalidParamsException {
    final File outFileDir = FileUtils.createTempDir("outputFilter", "test");
    try {
      final SequencesReader r = getReaderDNA();
      final ReaderParams rp = new MockReaderParams(r, SequenceMode.PROTEIN);
      final ISequenceParams sp = new MockSequenceParams(rp, 0, r.numberSequences());
      final NgsFilterParams filterParams = NgsFilterParams.builder().outputFilter(OutputFilter.PROTEIN_ALL_HITS).create();
      final NgsOutputParams ngsop = NgsOutputParams.builder().progress(false).outputDir(outFileDir).filterParams(filterParams).create();
      final NgsParams ngsp = NgsParams.builder().buildFirstParams(sp).searchParams(sp).outputParams(ngsop).proteinScoringMatrix(new ProteinScoringMatrix()).numberThreads(1).create();
      final OutputProcessor prot = ngsop.outFilter().makeProcessor(ngsp, null);
      prot.close();
      assertNotNull(prot);
      assertTrue(prot instanceof ProteinOutputProcessor);
    } finally {
      assertTrue(!outFileDir.exists() || FileHelper.deleteAll(outFileDir));
    }
  }

  public void testSingle1() throws IOException, InvalidParamsException {
    final File tmp = FileUtils.createTempDir("unfilteredope", "filter");
    try {
      final File input = new File(tmp, "1");
      ReaderTestUtils.getReaderDNA(">test\nacgtacgt", input, null).close();
      final File input2 = new File(tmp, "2");
      ReaderTestUtils.getReaderProtein(">test\nacgtacgt", input2).close();

      final NgsParams params = NgsParams.builder().outputParams(NgsOutputParams.builder().outputDir(tmp).create())
          .buildFirstParams(SequenceParams.builder().directory(input).useMemReader(true).create())
          .searchParams(SequenceParams.builder().loadNames(true).directory(input2).useMemReader(true).create())
          .proteinScoringMatrix(new ProteinScoringMatrix())
          .create();
      final OutputProcessor p = OutputFilter.SAM_SINGLE_END.makeProcessor(params, null);
      p.threadClone(HashingRegion.NONE).threadFinish();
      assertTrue(p instanceof SamSingleEndOutputProcessor);
      p.finish();
      p.close();
      assertTrue(new File(tmp, "alignments.sam").exists());
    } finally {
      assertTrue(FileHelper.deleteAll(tmp));
    }
  }

  public void testUnfiltered1() throws IOException, InvalidParamsException {
    final File tmp = FileUtils.createTempDir("unfilteredope", "filter");
    try {
      final File input = new File(tmp, "1");
      ReaderTestUtils.getReaderDNA(">test\nacgtacgt", input, null).close();
      final File input2 = new File(tmp, "2");
      ReaderTestUtils.getReaderProtein(">test\nacgtacgt", input2).close();

      final NgsParams params = NgsParams.builder().outputParams(NgsOutputParams.builder().outputDir(tmp).create())
          .buildFirstParams(SequenceParams.builder().directory(input).useMemReader(true).create())
          .searchParams(SequenceParams.builder().loadNames(true).directory(input2).useMemReader(true).create())
          .proteinScoringMatrix(new ProteinScoringMatrix())
          .create();
      final OutputProcessor p = OutputFilter.SAM_UNFILTERED.makeProcessor(params, null);
      p.threadClone(HashingRegion.NONE).threadFinish();
      assertTrue(p instanceof UnfilteredSingleEndOutputProcessor);
      p.finish();
      p.close();
      assertTrue(new File(tmp, "alignments.sam").exists());
    } finally {
      assertTrue(FileHelper.deleteAll(tmp));
    }
  }

  public void testUnfiltered2() throws IOException, InvalidParamsException {
    final File tmp = FileUtils.createTempDir("unfilteredse", "filter");
    try {
      final File input = new File(tmp, "1");
      ReaderTestUtils.getReaderDNA(">test\nacgtacgt", input, null).close();
      final File input2 = new File(tmp, "2");
      ReaderTestUtils.getReaderProtein(">test\nacgtacgt", input2).close();
      final File template = new File(tmp, "3");
      ReaderTestUtils.getReaderProtein(">test\nacgtacgt", template).close();

      try (NgsParams params = NgsParams.builder().outputParams(NgsOutputParams.builder().outputDir(tmp).create())
        .buildFirstParams(SequenceParams.builder().directory(input).useMemReader(true).create())
        .buildSecondParams(SequenceParams.builder().directory(input2).useMemReader(true).create())
        .maxFragmentLength(500).minFragmentLength(100)
        .searchParams(SequenceParams.builder().loadNames(true).directory(input2).useMemReader(true).create())
        .proteinScoringMatrix(new ProteinScoringMatrix())
        .create()) {
        try (OutputProcessor p = OutputFilter.SAM_UNFILTERED.makeProcessor(params, null)) {
          p.threadClone(HashingRegion.NONE).threadFinish();
          assertTrue(p instanceof UnfilteredPairedEndOutputProcessor);
          p.finish();
        }
      }
      assertTrue(new File(tmp, "alignments.sam").exists());
    } finally {
      assertTrue(FileHelper.deleteAll(tmp));
    }
  }


  public void testHashEquals() {
    TestUtils.equalsHashTest(new OutputFilter[][] {{OutputFilter.NONE}, {OutputFilter.PAIRED_END}, {OutputFilter.TOPN_PAIRED_END}});

  }

}
