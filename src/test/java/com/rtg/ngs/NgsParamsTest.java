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

import static com.rtg.util.StringUtils.LS;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;

import com.rtg.alignment.AlignerMode;
import com.rtg.launcher.ISequenceParams;
import com.rtg.launcher.SequenceParams;
import com.rtg.mode.ProgramMode;
import com.rtg.mode.SequenceMode;
import com.rtg.position.output.PositionParams;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.util.InvalidParamsException;
import com.rtg.util.MaxShiftFactor;
import com.rtg.util.StringUtils;
import com.rtg.util.TestUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.FileUtils;
import com.rtg.util.test.FileHelper;
import com.rtg.util.test.params.TestParams;

import junit.framework.TestCase;

/**
 */
public class NgsParamsTest extends TestCase {

  protected File mDir;

  @Override
  public void setUp() throws IOException {
    Diagnostic.setLogStream();
    mDir = FileHelper.createTempDirectory();
  }

  @Override
  public void tearDown() {
    assertTrue(FileHelper.deleteAll(mDir));
    mDir = null;
  }

  NgsParams getParams(final int threads, final SequenceParams build, final SequenceParams search, final NgsOutputParams outParams) {
    return NgsParams.builder().buildFirstParams(build)
        //                              .expectedInsertSize(7)
        .maxFragmentLength(8)
        .minFragmentLength(5)
        .numberThreads(threads)
        .outputParams(outParams)
        .maskParams("Bar")
        .searchParams(search)
        .stepSize(10)
        .compressHashes(true)
        .legacyCigars(true)
        .minHits(86)
        .intsetWindow(37)
        .enableProteinReadCache(true)
        .threadMultiplier(45)
        .readFreqThreshold(230)
        .mapXMinLength(23487)
        .parallelUnmatedProcessing(true)
        .substitutionPenalty(2)
        .gapOpenPenalty(3)
        .gapExtendPenalty(1)
        .alignerMode(AlignerMode.TABLE)
        .create();
  }

  public void testCloneBuider() throws Exception {
    final File subjectDir = ReaderTestUtils.getDNADir(mDir);
    final SequenceParams subjectaa = SequenceParams.builder().directory(subjectDir).mode(SequenceMode.UNIDIRECTIONAL).create();
    final File queryDir = ReaderTestUtils.getDNADir(mDir);
    final SequenceParams querya = SequenceParams.builder().directory(queryDir).create();
    final NgsFilterParams filterParamsa = NgsFilterParams.builder().outputFilter(OutputFilter.NONE).create();
    final NgsOutputParams outputa = NgsOutputParams.builder().progress(false).outputDir(new File("Foo")).filterParams(filterParamsa).create();
    final NgsParams params = getParams(3, subjectaa, querya, outputa);
    final NgsParams clone = params.cloneBuilder().create();

    //    assertEquals(params.expectedInsertSize(), clone.expectedInsertSize());
    assertEquals(params.buildFirstParams(), clone.buildFirstParams());
    assertEquals(params.listeners(), clone.listeners());
    assertEquals(params.maxFragmentLength(), clone.maxFragmentLength());
    assertEquals(params.minFragmentLength(), clone.minFragmentLength());
    assertEquals(params.name(), clone.name());
    assertEquals(params.numberThreads(), clone.numberThreads());
    assertEquals(params.outputParams(), clone.outputParams());
    assertEquals(params.proteinScoringMatrix(), clone.proteinScoringMatrix());
    assertEquals(params.searchParams(), clone.searchParams());
    assertEquals(params.buildSecondParams(), clone.buildSecondParams());
    assertEquals(params.stepSize(), clone.stepSize());
    assertEquals(params.useLongReadMapping(), clone.useLongReadMapping());
    assertEquals(params.compressHashes(), clone.compressHashes());

    assertEquals(params.legacyCigars(), clone.legacyCigars());
    assertEquals(params.minHits(), clone.minHits());
    assertEquals(params.intSetWindow(), clone.intSetWindow());
    assertEquals(params.enableProteinReadCache(), clone.enableProteinReadCache());
    assertEquals(params.threadMultiplier(), clone.threadMultiplier());
    assertEquals(params.readFreqThreshold(), clone.readFreqThreshold());
    assertEquals(params.mapXMinReadLength(), clone.mapXMinReadLength());
    assertEquals(params.parallelUnmatedProcessing(), clone.parallelUnmatedProcessing());

    assertEquals(params.gapOpenPenalty(), clone.gapOpenPenalty());
    assertEquals(params.gapExtendPenalty(), clone.gapExtendPenalty());
    assertEquals(params.substitutionPenalty(), clone.substitutionPenalty());
    assertEquals(params.alignerBandWidthFactor(), clone.alignerBandWidthFactor());
    assertEquals(params.alignerMode(), clone.alignerMode());
  }

  public void testEquals() throws Exception {
    final File subjectDir = ReaderTestUtils.getDNADir(mDir);
    final File subjectDir2 = ReaderTestUtils.getDNADir(mDir);
    final SequenceParams subjectaa = SequenceParams.builder().directory(subjectDir).mode(SequenceMode.UNIDIRECTIONAL).create();
    final SequenceParams subjectab = SequenceParams.builder().directory(subjectDir).mode(SequenceMode.UNIDIRECTIONAL).create();

    final SequenceParams subjectb = SequenceParams.builder().directory(subjectDir).mode(SequenceMode.UNIDIRECTIONAL).create();
    final SequenceParams subjectc = SequenceParams.builder().directory(subjectDir2).mode(SequenceMode.UNIDIRECTIONAL).create();

    final File queryDir = ReaderTestUtils.getDNADir(mDir);
    final SequenceParams querya = SequenceParams.builder().directory(queryDir).create();

    final File queryDirb = ReaderTestUtils.getDNADir(mDir);
    final SequenceParams queryba = SequenceParams.builder().directory(queryDirb).create();
    final File queryDirc = ReaderTestUtils.getDNADir(mDir);
    final SequenceParams querybb = SequenceParams.builder().directory(queryDirc).create();

    final NgsFilterParams filterParamsa = NgsFilterParams.builder().outputFilter(OutputFilter.NONE).create();
    final NgsOutputParams outputa = NgsOutputParams.builder().progress(false).outputDir(new File("Foo")).filterParams(filterParamsa).create();
    final NgsFilterParams filterParamsb = NgsFilterParams.builder().outputFilter(OutputFilter.PAIRED_END).create();
    final NgsOutputParams outputb = NgsOutputParams.builder().progress(false).outputDir(new File("Foo")).filterParams(filterParamsb).create();

    final NgsParams a1 = getParams(1, subjectaa, querya, outputa);
    final NgsParams a2 = getParams(1, subjectaa, querya, outputa);
    final NgsParams b = getParams(1, subjectaa, queryba, outputa);
    final NgsParams ba = getParams(1, subjectc, queryba, outputa);
    final NgsParams c = getParams(1, subjectb, queryba, outputb);
    final NgsParams d = getParams(1, subjectab, querybb, outputa);
    final NgsParams e = getParams(2, subjectab, querybb, outputa);
    TestUtils.equalsHashTest(new NgsParams[][] {{a1, a2}, {b}, {ba}, {c}, {d}, {e}});
    a1.close();
    a2.close();
    b.close();
    ba.close();
    c.close();
    d.close();
    e.close();
  }

  /** Subject sequence used for the calibration runs.  */
  public static final String SEQ_DNA_A1 = ""
      + ">x" + LS
      + "actg" + LS;

  /** Query sequence used for the calibration runs.  */
  public static final String SEQ_DNA_A2 = ""
      + ">u" + LS
      + "actgact" + LS
      + ">v" + LS
      + "antg" + LS;

  public void test() throws Exception {
    final File testFile = FileHelper.createTempFile();
    try {
      assertTrue(testFile.delete());
      final ProgramMode pm = ProgramMode.SLIMN;
      final File subjectDir = ReaderTestUtils.getDNASubDir(SEQ_DNA_A1, mDir);
      final SequenceParams subject = SequenceParams.builder().directory(subjectDir).mode(pm.subjectMode()).create();

      final NgsFilterParams filterParams = NgsFilterParams.builder().outputFilter(OutputFilter.NONE).create();
      final NgsOutputParams count = NgsOutputParams.builder().progress(false).outputDir(testFile).filterParams(filterParams).create();
      final File queryDir = ReaderTestUtils.getDNASubDir(SEQ_DNA_A2, mDir);
      final SequenceParams query = SequenceParams.builder().directory(queryDir).mode(pm.queryMode()).create();

      final NgsParams bsp = getParams(3, subject, query, count);
      bsp.integrity();
      assertEquals(new File(testFile, "bar"), bsp.outputParams().file("bar"));
      assertEquals(3, bsp.numberThreads());
      assertFalse(bsp.paired());
      assertEquals(Integer.valueOf(8), bsp.maxFragmentLength());
      assertEquals(Integer.valueOf(5), bsp.minFragmentLength());
      assertTrue(bsp.parallelUnmatedProcessing());
      assertNull(bsp.buildSecondParams());
      assertEquals(AlignerMode.TABLE, bsp.alignerMode());
      final File file = bsp.file("foobar");
      assertFalse(file == null);
      final String name = file.toString();
      assertTrue(name.endsWith(testFile.getPath() + StringUtils.FS + "foobar"));
      final OutputStream stream = bsp.unusedOutStream();
      try {
        assertFalse(stream == null);
        stream.close();
        assertFalse(bsp.listeners() == null);
        assertTrue(bsp.directory().toString().endsWith(testFile.getPath()));
        assertFalse(bsp.directory() == null);
        try {
          bsp.integrity();

          assertEquals(subject.toString(), bsp.buildFirstParams().toString());
          assertEquals(query.toString(), bsp.searchParams().toString());

          bsp.close();
          bsp.buildFirstParams().reader();
          bsp.close();
          bsp.searchParams().reader();
          bsp.close();
          bsp.buildFirstParams().reader();
          bsp.searchParams().reader();
        } finally {
          assertTrue(FileHelper.deleteAll(bsp.directory()));
        }
      } finally {
        bsp.close();
      }
    } finally {
      assertTrue(!testFile.exists() || FileHelper.deleteAll(testFile));
    }
  }


  public void testCompressOutput() throws Exception {
    final File testFile = FileHelper.createTempFile();
    try {
      assertTrue(testFile.delete());
      final ProgramMode pm = ProgramMode.SLIMN;
      final File subjectDir = ReaderTestUtils.getDNADir(mDir);
      final SequenceParams subject = SequenceParams.builder().directory(subjectDir).mode(pm.subjectMode()).create();

      final NgsFilterParams filterParams = NgsFilterParams.builder().outputFilter(OutputFilter.NONE).zip(true).create();
      final NgsOutputParams outputParams = NgsOutputParams.builder().progress(false).outputDir(testFile).filterParams(filterParams).create();
      final File queryDir = ReaderTestUtils.getDNADir(mDir);
      final SequenceParams query = SequenceParams.builder().directory(queryDir).mode(pm.queryMode()).create();

      try (NgsParams bsp = getParams(3, subject, query, outputParams)) {
        assertTrue(bsp.compressOutput());
      }
    } finally {
      assertTrue(!testFile.exists() || FileHelper.deleteAll(testFile));
    }
  }

  NgsParams makeParams(int w, int a, int b, int c) throws IOException, InvalidParamsException {
    final File sub = ReaderTestUtils.getDNADir(mDir); //FileUtils.createTempDir("sub", "ngs");
    final File que = ReaderTestUtils.getDNADir(mDir);
    final File res = FileUtils.createTempDir("res", "ngs");
    FileHelper.deleteAll(res);
    final NgsParamsBuilder builder = NgsParams.builder();
    builder.searchParams(SequenceParams.builder().directory(que).mode(SequenceMode.UNIDIRECTIONAL).useMemReader(true).create());
    builder.buildFirstParams(SequenceParams.builder().directory(sub).mode(SequenceMode.UNIDIRECTIONAL).useMemReader(true).create());
    builder.outputParams(NgsOutputParams.builder().outputDir(res)
        .create());
    builder.maskParams(new NgsMaskParamsGeneral(w, a, b, c));
    builder.numberThreads(1);
    builder.stepSize(w);
    return builder.create();
  }

  void closeParams(final NgsParams params) throws IOException {
    params.close();
    FileHelper.deleteAll(params.searchParams().directory());
    FileHelper.deleteAll(params.buildFirstParams().directory());
    final ISequenceParams buildSecond = params.buildSecondParams();
    if (buildSecond != null) {
      FileHelper.deleteAll(buildSecond.directory());
    }
  }

  public void testToPositionParams1() throws IOException, InvalidParamsException {
    final NgsParams params = makeParams(2, 1, 2, 1);
    final PositionParams pparams = params.toPositionParams();

    assertEquals(2, pparams.output().distribution().maxIndel());

    //System.err.println(params.toString());
    final String pps = pparams.toString();
    //System.err.println(pps);
    TestUtils.containsAll(pps, "distribution={maxGap=4}}",
      "hash bits=4 initial pointer bits=2 value bits=1",
      "window=2 step=2",
      "window=2 step=1");
    closeParams(params);
  }

  //maxGap = 0 (computed from r, w, a, b, c)
  public void testToPositionParams2() throws IOException, InvalidParamsException {
    final NgsParams params = makeParams(4, 0, 0, 1);
    final PositionParams pparams = params.toPositionParams();
    //System.err.println(params.toString());
    final String pps = pparams.toString();
    //System.err.println(pps);
    TestUtils.containsAll(pps, "distribution={maxGap=1}}",
      "hash bits=8 initial pointer bits=2 value bits=0",
      "window=4 step=4",
      "window=4 step=1");
    closeParams(params);
  }

  public void testOmnes() {
    new TestParams(NgsParams.class, NgsParamsBuilder.class).check();
  }

  public void testAlignerBandWidthFactor() throws Exception {
    final File testFile = FileHelper.createTempFile();
    try {
      assertTrue(testFile.delete());
      final ProgramMode pm = ProgramMode.SLIMN;
      final File subjectDir = ReaderTestUtils.getDNADir(mDir);
      final SequenceParams subject = SequenceParams.builder().directory(subjectDir).mode(pm.subjectMode()).create();

      final NgsFilterParams filterParams = NgsFilterParams.builder().outputFilter(OutputFilter.NONE).create();
      final NgsOutputParams count = NgsOutputParams.builder().progress(false).outputDir(testFile).filterParams(filterParams).create();
      final File queryDir = ReaderTestUtils.getDNADir(mDir);
      final SequenceParams query = SequenceParams.builder().directory(queryDir).mode(pm.queryMode()).create();

      final NgsParams params = getParams(3, subject, query, count);

      assertEquals(2, params.substitutionPenalty());
      assertEquals(3, params.gapOpenPenalty());
      assertEquals(1, params.gapExtendPenalty());
      assertNotNull(params.alignerBandWidthFactor());
      assertEquals(0.5, params.alignerBandWidthFactor().getFactor(), 0.0001);
      
      final NgsParams params2 = params.cloneBuilder().alignerBandWidthFactor(new MaxShiftFactor(0.4)).create();

      assertEquals(2, params2.substitutionPenalty());
      assertEquals(3, params2.gapOpenPenalty());
      assertEquals(1, params2.gapExtendPenalty());
      assertNotNull(params2.alignerBandWidthFactor());
      assertEquals(0.4, params2.alignerBandWidthFactor().getFactor(), 0.0001);

    } finally {
      assertTrue(!testFile.exists() || FileHelper.deleteAll(testFile));
    }

  }
}
