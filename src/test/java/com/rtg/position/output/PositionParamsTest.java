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
package com.rtg.position.output;

import static com.rtg.util.StringUtils.LS;

import java.io.File;
import java.io.IOException;

import com.rtg.index.IndexUtils;
import com.rtg.index.params.CreateParams;
import com.rtg.launcher.BuildParams;
import com.rtg.launcher.HashingRegion;
import com.rtg.launcher.ISequenceParams;
import com.rtg.launcher.MockReaderParams;
import com.rtg.launcher.MockSequenceParams;
import com.rtg.launcher.ReaderParams;
import com.rtg.launcher.SequenceParams;
import com.rtg.mode.ProgramMode;
import com.rtg.position.PositionUtils;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.util.TestUtils;
import com.rtg.util.Utils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.test.FileHelper;

import junit.framework.TestCase;

/**
 */
public class PositionParamsTest extends TestCase {

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

  /**
   * Compute total number of bytes used by <code>BuildSearch</code>
   * @param buildSearchParams parameters to be used to create <code>BuildSearch</code>
   * @return total number of bytes used by <code>BuildSearch</code>
   */
  static long bytes(final PositionParams buildSearchParams) {
    return IndexUtils.bytes(buildSearchParams.indexParams());
  }

  PositionParams getParams(final ProgramMode mode, final Integer threshold, final BuildParams build, final BuildParams search, final PositionOutputParams outParams, final int numberThreads, final boolean largeSort) throws IOException {
    return PositionParams.builder().mode(mode).hashCountThreshold(threshold).buildParams(build).searchParams(search).indexParams(CreateParams.fromBuildParams(build).create()).outputParams(outParams).numberThreads(numberThreads).create();
  }

  public void testLog10() {
    assertEquals(1, PositionParams.log10(1));
    assertEquals(1, PositionParams.log10(9));
    assertEquals(1, PositionParams.log10(10));
    assertEquals(2, PositionParams.log10(11));
    assertEquals(2, PositionParams.log10(99));
    assertEquals(2, PositionParams.log10(100));
    assertEquals(3, PositionParams.log10(101));
    assertEquals(4, PositionParams.log10(1001));
    assertEquals(10, PositionParams.log10(Integer.MAX_VALUE));
    try {
      PositionParams.log10(0);
      fail();
    } catch (final RuntimeException e) {
      //expected
    }
  }

  public void testZeroFormat() {
    assertEquals("", PositionParams.zeroFormat(0, 1));
    assertEquals("0", PositionParams.zeroFormat(0, 2));
    assertEquals("9", PositionParams.zeroFormat(9, 10));
    assertEquals("09", PositionParams.zeroFormat(9, 11));
    assertEquals("99", PositionParams.zeroFormat(99, 100));
    assertEquals("099", PositionParams.zeroFormat(99, 101));
    assertEquals("2147483646", PositionParams.zeroFormat(Integer.MAX_VALUE - 1, Integer.MAX_VALUE));
  }

  public void testHash() throws IOException {
    final ProgramMode pma = ProgramMode.SLIMN;
    final ReaderParams rp = new MockReaderParams(0, 0, pma.subjectMode().codeType());
    final ISequenceParams subjectaa = new MockSequenceParams(rp, pma.subjectMode(), 0, 0);
    final BuildParams buildaa = BuildParams.builder().windowSize(4).stepSize(1).sequences(subjectaa).create();

    final File hitDir = new File("");

    final ReaderParams qp = new MockReaderParams(0, 0, pma.queryMode().codeType());
    final ISequenceParams querya = new MockSequenceParams(qp, pma.queryMode(), 0, 0);
    final BuildParams queriesa = BuildParams.builder().windowSize(4).stepSize(1).sequences(querya).create();

    final PositionOutputParams outParams = new PositionOutputParams(hitDir, OutputFormatType.SEGMENT, null, null, false, 1);
    final PositionParams a1 = PositionParams.builder().mode(pma).buildParams(buildaa).searchParams(queriesa).outputParams(outParams).create();
    final PositionParams a2 = PositionParams.builder().mode(pma).hashCountThreshold(100).buildParams(buildaa).searchParams(queriesa).outputParams(outParams).create();
    assertTrue(a1.hashCode() != a2.hashCode());
  }

  public void testEquals() throws IOException {
    final ProgramMode pma = ProgramMode.SLIMN;
    final ProgramMode pmb = ProgramMode.TSLIMX;

    final File subjectDir = ReaderTestUtils.getDNADir(mDir);
    final SequenceParams subjectaa = SequenceParams.builder().directory(subjectDir).mode(pma.subjectMode()).create();
    final BuildParams buildaa = BuildParams.builder().windowSize(4).stepSize(1).sequences(subjectaa).create();
    final SequenceParams subjectab = SequenceParams.builder().directory(subjectDir).mode(pmb.subjectMode()).create();
    final BuildParams buildab = BuildParams.builder().windowSize(4).stepSize(2).sequences(subjectab).create();

    final SequenceParams subjectb = SequenceParams.builder().directory(subjectDir).mode(pmb.subjectMode()).create();
    final BuildParams buildb = BuildParams.builder().windowSize(4).stepSize(1).sequences(subjectb).create();

    final File hitDir = FileHelper.createTempDirectory(mDir);
    final File queryDir = ReaderTestUtils.getDNASubDir(SEQ_DNA_A2, mDir);
    final SequenceParams querya = SequenceParams.builder().directory(queryDir).mode(pma.queryMode()).create();
    final BuildParams queriesa = BuildParams.builder().windowSize(4).stepSize(1).sequences(querya).create();

    final SequenceParams queryba = SequenceParams.builder().directory(queryDir).mode(pmb.queryMode()).create();
    final BuildParams queriesba = BuildParams.builder().windowSize(4).stepSize(1).sequences(queryba).create();
    final SequenceParams querybb = SequenceParams.builder().directory(queryDir).mode(pmb.queryMode()).create();
    final BuildParams queriesbb = BuildParams.builder().windowSize(5).stepSize(1).sequences(querybb).create();

    final PositionOutputParams outParamsa = new PositionOutputParams(hitDir, OutputFormatType.SEGMENT, null, null, false, 1);
    final PositionOutputParams outParamsb = new PositionOutputParams(hitDir, OutputFormatType.SEGMENT, null, null, true, 1);
    final PositionParams a1 = PositionParams.builder().mode(pma).buildParams(buildaa).searchParams(queriesa).outputParams(outParamsa).create();
    final PositionParams a2 = PositionParams.builder().mode(pma).buildParams(buildaa).searchParams(queriesa).outputParams(outParamsa).create();
    final PositionParams b = PositionParams.builder().mode(pmb).buildParams(buildb).searchParams(queriesba).outputParams(outParamsa).create();
    final PositionParams c = PositionParams.builder().mode(pmb).buildParams(buildab).searchParams(queriesba).outputParams(outParamsa).create();
    final PositionParams d = PositionParams.builder().mode(pmb).buildParams(buildab).searchParams(queriesbb).outputParams(outParamsa).create();
    final PositionParams e = PositionParams.builder().mode(pmb).hashCountThreshold(1).buildParams(buildab).searchParams(queriesbb).outputParams(outParamsa).create();
    final PositionParams f = PositionParams.builder().mode(pmb).hashCountThreshold(2).buildParams(buildab).searchParams(queriesbb).outputParams(outParamsa).create();
    final PositionParams g = PositionParams.builder().mode(pmb).hashCountThreshold(2).buildParams(buildab).searchParams(queriesbb).outputParams(outParamsa).progress(true).create();
    final PositionParams j = PositionParams.builder().mode(pmb).hashCountThreshold(2).buildParams(buildab).searchParams(queriesbb).outputParams(outParamsb).progress(true).numberThreads(2/*threads*/).create();
    TestUtils.equalsHashTest(new PositionParams[][] {{a1, a2}, {b}, {c}, {d}, {e}, {f}, {g}, {j}});
    a1.close();
    a2.close();
    b.close();
    c.close();
    d.close();
    e.close();
    f.close();
    g.close();
    j.close();
  }

  /** Query sequence used for the calibration runs.  */
  public static final String SEQ_DNA_A2 = ""
    + ">u" + LS
    + "actgact" + LS
    + ">v" + LS
    + "antg" + LS;

  /** Query sequence used for the calibration runs.  */
  public static final String SEQ_DNA_A3 = ""
    + ">x1" + LS
    + "actgact" + LS
    + ">x2" + LS
    + "actgact" + LS
    + ">x3" + LS
    + "actgact" + LS
    + ">x4" + LS
    + "actgact" + LS
    ;

  public void test() throws Exception {
    final ProgramMode pm = ProgramMode.SLIMN;
    final File subjectDir = ReaderTestUtils.getDNADir(mDir);
    final SequenceParams subject = SequenceParams.builder().directory(subjectDir).mode(pm.subjectMode()).create();
    final BuildParams build = BuildParams.builder().windowSize(4).stepSize(1).sequences(subject).create();

    final File hitDir = FileHelper.createTempDirectory(mDir);
    final File queryDir = ReaderTestUtils.getDNASubDir(SEQ_DNA_A2, mDir);
    final SequenceParams query = SequenceParams.builder().directory(queryDir).mode(pm.queryMode()).create();
    final BuildParams queries = BuildParams.builder().windowSize(4).stepSize(1).sequences(query).create();
    final PositionOutputParams outputParams = new PositionOutputParams(hitDir, OutputFormatType.SEGMENT, null, null, false, 1);
    final PositionParams wp = getParams(pm, 1000, build, queries, outputParams, 13/*thread*/, false/*largeSort*/);
    assertEquals(new File(hitDir.getPath(), "log").getPath(), wp.logFile());
    try {
      wp.integrity();

      assertEquals(pm, wp.mode());
      assertEquals(Integer.valueOf(1000), wp.hashCountThreshold());
      assertFalse(wp.progress());
      assertEquals(build.toString(), wp.build().toString());
      assertEquals(queries.toString(), wp.search().toString());
      assertEquals(7, wp.bufferLength());
      assertEquals(13, wp.numberThreads());
      assertEquals(15, bytes(wp));
      assertEquals(""
          + "PositionParams mode=SLIMN threshold=1000 progress=" + Boolean.FALSE.toString() + " number threads=13" + LS
          + ".. output={outputDir=" + wp.output().directory() + " format=SEGMENT zip=" + Boolean.FALSE.toString()  + " score threshold=" + Utils.realFormat(null) + " topN=1" + " distribution={" + Utils.realFormat(null) + "}}" + LS
          + ".. index={ size=1 hash bits=8 initial pointer bits=2 value bits=20}" + LS
          + ".. build={ seq={SequenceParams mode=UNIDIRECTIONAL region=[(0:-1), (1:-1)] directory="
          + build.directory().toString()
          + "}  window=4 step=1}" + LS
          + ".. search={ seq={SequenceParams mode=BIDIRECTIONAL region=[(0:-1), (2:-1)] directory="
          + queries.directory()
          + "}  window=4 step=1}"  + LS
          , wp.toString()
      );

      assertEquals(""
          + "\tMemory\tShared_buffer\t7" + LS
          + "\tMemory\tHash\t1" + LS
          + "\tMemory\tValue\t4" + LS
          + "\tMemory\tInitial_position\t6" + LS
          + "\tMemory\tBit_vector\t4" + LS
          , PositionUtils.memToString(wp)
      );
      wp.close();
      wp.build().sequences().reader();
      wp.close();
      wp.search().sequences().reader();
      wp.close();
      wp.build().sequences().reader();
      wp.search().sequences().reader();
    } finally {
      wp.close();
    }
  }

  public void testSearchSubSequence() throws Exception {
    final ProgramMode pm = ProgramMode.SLIMN;
    final File subjectDir = ReaderTestUtils.getDNADir(mDir);
    final SequenceParams subject = SequenceParams.builder().directory(subjectDir).mode(pm.subjectMode()).create();
    final BuildParams build = BuildParams.builder().windowSize(4).stepSize(1).sequences(subject).create();

    final File hitDir = FileHelper.createTempDirectory(mDir);
    final File queryDir = ReaderTestUtils.getDNASubDir(SEQ_DNA_A3, mDir);
    final SequenceParams query = SequenceParams.builder().directory(queryDir).mode(pm.queryMode()).create();
    final BuildParams queries = BuildParams.builder().windowSize(4).stepSize(1).sequences(query).create();
    final PositionOutputParams outputParams = new PositionOutputParams(hitDir, OutputFormatType.SEGMENT, null, null, false, 1);

    final PositionParams pp1 = getParams(pm, 1, build, queries, outputParams, 1/*thread*/, false/*largeSort*/);
    pp1.integrity();
    final PositionParams pp2 = pp1.subSearch(new HashingRegion(1, 2));
    final String first = "PositionParams mode=SLIMN threshold=1 progress=" + Boolean.FALSE.toString() + " number threads=1" + LS
    + ".. output={outputDir=" + pp1.outputDir() + " format=SEGMENT zip=" + Boolean.FALSE.toString() + " score threshold=" + Utils.realFormat(null) + " topN=1 distribution={" + Utils.realFormat(null) + "}}" + LS
    + ".. index={ size=1 hash bits=8 initial pointer bits=2 value bits=20}" + LS
    + ".. build={ seq={SequenceParams mode=UNIDIRECTIONAL region=[(0:-1), (1:-1)] directory="
    + build.directory().toString()
    + "}  window=4 step=1}" + LS;
    assertEquals(""
        + first
        + ".. search={ seq={SequenceParams mode=BIDIRECTIONAL region=[(0:-1), (4:-1)] directory=" + queries.directory() + "}  window=4 step=1}" + LS
        , pp1.toString()
    );
    assertEquals(""
        + first
        + ".. search={ seq={SequenceParams mode=BIDIRECTIONAL region=[(1:-1), (2:-1)] directory=" + queries.directory() + "}  window=4 step=1}" + LS
        , pp2.toString()
    );
    pp1.close();
    pp2.close();
  }

  public void testBuildSubSequence() throws Exception {
    final ProgramMode pm = ProgramMode.SLIMN;
    final File subjectDir = ReaderTestUtils.getDNADir(mDir);
    final SequenceParams subject = SequenceParams.builder().directory(subjectDir).mode(pm.subjectMode()).create();
    final BuildParams build = BuildParams.builder().windowSize(4).stepSize(1).sequences(subject).create();

    final File hitDir = FileHelper.createTempDirectory(mDir);
    final File queryDir = ReaderTestUtils.getDNASubDir(SEQ_DNA_A3, mDir);
    final SequenceParams query = SequenceParams.builder().directory(queryDir).mode(pm.queryMode()).create();
    final BuildParams queries = BuildParams.builder().windowSize(4).stepSize(1).sequences(query).create();
    final PositionOutputParams outputParams = new PositionOutputParams(hitDir, OutputFormatType.SEGMENT, null, null, false, 1);

    final PositionParams pp1 = getParams(pm, 1, build, queries, outputParams, 1/*thread*/, false/*largeSort*/);
    pp1.integrity();
    final PositionParams pp2 = pp1.subBuild(new HashingRegion(0, 1));
    final String indexStr = ".. index={ size=1 hash bits=8 initial pointer bits=2 value bits=20}" + LS;
    final String buildStr = ".. build={ seq={SequenceParams mode=UNIDIRECTIONAL region=[(0:-1), (1:-1)] directory=";
    final String secondStr = build.directory().toString() + "}  window=4 step=1}" + LS;
    final String first = "PositionParams mode=SLIMN threshold=1 progress=" + Boolean.FALSE.toString() + " number threads=1" + LS + ".. output={outputDir=" + pp1.outputDir() + " format=SEGMENT zip=" + Boolean.FALSE.toString() + " score threshold=" + Utils.realFormat(null) + " topN=1 distribution={" + Utils.realFormat(null) + "}}" + LS;
    final String searchStr = ".. search={ seq={SequenceParams mode=BIDIRECTIONAL region=[(0:-1), (4:-1)] directory=" + queries.directory() + "}  window=4 step=1}" + LS;
    assertEquals(first + indexStr + buildStr + secondStr + searchStr, pp1.toString());
    final String buildStr2 = ".. build={ seq={SequenceParams mode=UNIDIRECTIONAL region=[(0:-1), (1:-1)] directory=";
    assertEquals(first + indexStr + buildStr2 + secondStr + searchStr, pp2.toString());
    pp1.close();
    pp2.close();
  }

  public void test1() throws Exception {
    final ProgramMode pm = ProgramMode.SLIMN;
    final File subjectDir = ReaderTestUtils.getDNADir(mDir);
    final SequenceParams subject = SequenceParams.builder().directory(subjectDir).mode(pm.subjectMode()).create();
    final BuildParams build = BuildParams.builder().windowSize(4).stepSize(1).sequences(subject).create();

    final File hitDir = FileHelper.createTempDirectory(mDir);
    final File queryDir = ReaderTestUtils.getDNASubDir(SEQ_DNA_A2, mDir);
    final SequenceParams query = SequenceParams.builder().directory(queryDir).mode(pm.queryMode()).create();
    final BuildParams queries = BuildParams.builder().windowSize(4).stepSize(1).sequences(query).create();
    final PositionOutputParams outputParams = new PositionOutputParams(hitDir, OutputFormatType.SEGMENT, null, null, false, 1);

    final PositionParams bsp = getParams(pm, 1, build, queries, outputParams, 1/*thread*/, true/*largeSort*/);
    assertEquals(new File(hitDir.getPath(), "log").getPath(), bsp.logFile());
    try {
      bsp.integrity();

      assertEquals(pm, bsp.mode());
      assertEquals(Integer.valueOf(1), bsp.hashCountThreshold());
      assertEquals(build.toString(), bsp.build().toString());
      assertEquals(queries.toString(), bsp.search().toString());
      assertEquals(7, bsp.bufferLength());
      assertEquals(1, bsp.numberThreads());
      assertEquals(15, bytes(bsp));
      assertEquals(""
          + "PositionParams mode=SLIMN threshold=1 progress=" + Boolean.FALSE.toString() + " number threads=1" + LS
          + ".. output={outputDir=" + bsp.outputDir() + " format=SEGMENT zip=" + Boolean.FALSE.toString() + " score threshold=" + Utils.realFormat(null) + " topN=1 distribution={" + Utils.realFormat(null) + "}}" + LS
          + ".. index={ size=1 hash bits=8 initial pointer bits=2 value bits=20}" + LS
          + ".. build={ seq={SequenceParams mode=UNIDIRECTIONAL region=[(0:-1), (1:-1)] directory="
          + build.directory().toString()
          + "}  window=4 step=1}" + LS
          + ".. search={ seq={SequenceParams mode=BIDIRECTIONAL region=[(0:-1), (2:-1)] directory="
          + queries.directory()
          + "}  window=4 step=1}"  + LS
          , bsp.toString()
      );

      assertEquals(""
          + "\tMemory\tShared_buffer\t7" + LS
          + "\tMemory\tHash\t1" + LS
          + "\tMemory\tValue\t4" + LS
          + "\tMemory\tInitial_position\t6" + LS
          + "\tMemory\tBit_vector\t4" + LS
          , PositionUtils.memToString(bsp)
      );
      bsp.close();
      bsp.build().sequences().reader();
      bsp.close();
      bsp.search().sequences().reader();
      bsp.close();
      bsp.build().sequences().reader();
      bsp.search().sequences().reader();
    } finally {
      bsp.close();
    }
  }

}

