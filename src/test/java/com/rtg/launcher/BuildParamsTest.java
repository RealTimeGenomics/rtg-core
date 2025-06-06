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
package com.rtg.launcher;

import java.io.File;
import java.io.IOException;

import com.rtg.mode.ProgramMode;
import com.rtg.mode.SequenceMode;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.util.StringUtils;
import com.rtg.util.TestUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.intervals.LongRange;
import com.rtg.util.io.FileUtils;
import com.rtg.util.test.FileHelper;

import junit.framework.TestCase;

/**
 */
public class BuildParamsTest extends TestCase {

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

  private static final String SEQ_DNA_A = ""
    + ">x" + StringUtils.LS
    + "actg" + StringUtils.LS;

  private static final String SEQ_DNA_B = ""
    + ">x1" + StringUtils.LS
    + "actg" + StringUtils.LS
    + ">x2" + StringUtils.LS
    + "actg" + StringUtils.LS
    + ">x3" + StringUtils.LS
    + "actg" + StringUtils.LS
    ;

  BuildParams getParams(final SequenceMode mode, final File subjectsDir, final int windowSize, final int stepSize) throws Exception {
    final SequenceParams seq = SequenceParams.builder().directory(subjectsDir).mode(mode).create();
    return BuildParams.builder().windowSize(windowSize).stepSize(stepSize).sequences(seq).create();
  }

  BuildParams getParams(final long size, final int windowSize, final int stepSize) throws Exception {
    final File file = ReaderTestUtils.getDNADir(mDir);
    final SequenceParams sp = SequenceParams.builder().directory(file).create();
    try {
      return BuildParams.builder().windowSize(windowSize).stepSize(stepSize).sequences(sp).create();
    } catch (final RuntimeException re) {
      sp.close();
      throw re;
    }
  }

  public void testSubSequence() throws IOException {
    final File filea = ReaderTestUtils.getDNASubDir(SEQ_DNA_B, mDir);
    final SequenceParams spa = SequenceParams.builder().directory(filea).create();
    spa.close();
    final BuildParams a = BuildParams.builder().windowSize(4).stepSize(2).sequences(spa).create();
    //System.err.println(a.toString());
    final BuildParams b = a.subSequence(new HashingRegion(1, 2));
    final String actual = b.toString();
    //System.err.println(actual);
    assertTrue(actual, actual.startsWith(" seq={SequenceParams mode=BIDIRECTIONAL region=[(1:-1), (2:-1)] directory="));
    assertTrue(actual, actual.endsWith("}  window=4 step=2"));
    a.close();
    b.close();
  }

  public void testEquals() throws IOException {
    final File filea = ReaderTestUtils.getDNADir(mDir);
    final SequenceParams spa = SequenceParams.builder().directory(filea).create();
    spa.close();
    final File fileb = ReaderTestUtils.getDNADir(mDir);
    final SequenceParams spb = SequenceParams.builder().directory(fileb).create();
    spb.close();
    final BuildParams a1 = BuildParams.builder().windowSize(4).stepSize(2).sequences(spa).create();
    final BuildParams a2 = BuildParams.builder().windowSize(4).stepSize(2).sequences(spa).create();
    final BuildParams d = BuildParams.builder().windowSize(5).stepSize(2).sequences(spa).create();
    final BuildParams e = BuildParams.builder().windowSize(4).stepSize(3).sequences(spa).create();
    TestUtils.equalsHashTest(new BuildParams[][] {{a1, a2}, {d}, {e}});
    a1.close();
    a2.close();
    d.close();
    e.close();
  }

  private static final String MEM_EXPECTED = ""
    + "\tMemory\tHash\t42" + StringUtils.LS
    + "\tMemory\tValue\t168" + StringUtils.LS
    + "\tMemory\tInitial_position\t34" + StringUtils.LS
    + "\tMemory\tBit_vector\t32" + StringUtils.LS;

  public void testa() throws Exception {
    final BuildParams ip = getParams(42L, 4, 3);
    assertEquals(4, ip.windowSize());
    assertEquals(3, ip.stepSize());
    final String str = ip.toString();
    assertTrue(str, str.endsWith(" window=4 step=3"));

    ip.sequences().reader();
    ip.close();
  }

  public void testHashBits0() throws Exception {
    final BuildParams ip = getParams(42L, 33, 3);
    assertEquals(33, ip.windowSize());
    assertEquals(3, ip.stepSize());
    final String str = ip.toString();
    assertTrue(str, str.endsWith(" window=33 step=3"));
    ip.close();
  }

  public void testHashBits1() throws Exception {
    final BuildParams ip = getParams(42L, 32, 3);
    assertEquals(32, ip.windowSize());
    assertEquals(3, ip.stepSize());
    final String str = ip.toString();
    assertTrue(str, str.endsWith(" window=32 step=3"));
    ip.close();
  }

  public void testHashBits2() throws Exception {
    final BuildParams ip = getParams(42L, 13, 3);
    assertEquals(13, ip.windowSize());
    assertEquals(3, ip.stepSize());
    final String str = ip.toString();
    assertTrue(str, str.endsWith(" window=13 step=3"));
    ip.close();
  }


  public void testBasic1() throws Exception {
    File file = null;
    try {
      file = ReaderTestUtils.getDNADir(mDir);
      try (BuildParams ip = getParams(ProgramMode.SLIMN.subjectMode(), file, 4, 3)) {
        assertEquals(1, ip.calculateSize());
        assertEquals(4, ip.windowSize());
        assertEquals(3, ip.stepSize());
        assertEquals(file, ip.directory());
        final String str = ip.toString();
        //System.err.println(str);
        assertTrue(str, str.startsWith(" seq={SequenceParams mode="));
        assertTrue(str, str.endsWith(" window=4 step=3"));
      }
    } finally {
      FileHelper.deleteAll(file);
    }
  }

  public void testBasic2() throws Exception {
    File file = null;
    try {
      file = ReaderTestUtils.getDNADir(mDir);
      try (BuildParams ip = getParams(ProgramMode.TSLIMX.subjectMode(), file, 4, 3)) {
        assertEquals(0, ip.calculateSize());
        assertEquals(4, ip.windowSize());
        assertEquals(3, ip.stepSize());
        assertEquals(file, ip.directory());
        final String str = ip.toString();
        //System.err.println(str);
        assertTrue(str, str.startsWith(" seq={SequenceParams mode="));
        assertTrue(str, str.endsWith(" window=4 step=3"));
      }
    } finally {
      FileHelper.deleteAll(file);
    }
  }

  public void testSize() {
    final long numSeqs = 10;
    assertEquals(0, BuildParams.size(0, numSeqs, 1, 1, 1, 1));
    assertEquals(1, BuildParams.size(1, numSeqs, 1, 1, 1, 1));
    assertEquals(0, BuildParams.size(0, numSeqs, 5, 2, 1, 1));
    assertEquals(0, BuildParams.size(1, numSeqs, 5, 2, 1, 1));
    assertEquals(0, BuildParams.size(2, numSeqs, 5, 2, 1, 1));
    assertEquals(0, BuildParams.size(3, numSeqs, 5, 2, 1, 1));
    assertEquals(0, BuildParams.size(4, numSeqs, 5, 2, 1, 1));
    assertEquals(1, BuildParams.size(5, numSeqs, 5, 2, 1, 1));
    assertEquals(1, BuildParams.size(6, numSeqs, 5, 2, 1, 1));
    assertEquals(2, BuildParams.size(7, numSeqs, 5, 2, 1, 1));
    assertEquals(2, BuildParams.size(8, numSeqs, 5, 2, 1, 1));
    assertEquals(4, BuildParams.size(8, numSeqs, 5, 2, 2, 1));
    assertEquals(0, BuildParams.size(8, numSeqs, 5, 2, 2, 3));
    assertEquals(4, BuildParams.size(24, numSeqs, 5, 2, 2, 3));

    assertEquals(24, BuildParams.size(1000, 10, 25, 75, 1, 1)); // Check adding extra in when step > window
  }

  public void testSizeBuild() throws IOException {
    assertEquals(0, sizeCheck(0, 1, 1, SequenceMode.BIDIRECTIONAL));
    assertEquals(1, sizeCheck(1, 1, 1, SequenceMode.UNIDIRECTIONAL));
    assertEquals(2, sizeCheck(1, 1, 1, SequenceMode.BIDIRECTIONAL));
    assertEquals(0, sizeCheck(0, 5, 2, SequenceMode.BIDIRECTIONAL));
    assertEquals(0, sizeCheck(1, 5, 2, SequenceMode.BIDIRECTIONAL));
    assertEquals(0, sizeCheck(2, 5, 2, SequenceMode.BIDIRECTIONAL));
    assertEquals(0, sizeCheck(3, 5, 2, SequenceMode.BIDIRECTIONAL));
    assertEquals(0, sizeCheck(4, 5, 2, SequenceMode.BIDIRECTIONAL));
    assertEquals(1, sizeCheck(5, 5, 2, SequenceMode.UNIDIRECTIONAL));
    assertEquals(1, sizeCheck(5, 5, 2, SequenceMode.PROTEIN));
    assertEquals(2, sizeCheck(5, 5, 2, SequenceMode.BIDIRECTIONAL));
    assertEquals(2, sizeCheck(6, 5, 2, SequenceMode.BIDIRECTIONAL));
    assertEquals(4, sizeCheck(7, 5, 2, SequenceMode.BIDIRECTIONAL));
    assertEquals(4, sizeCheck(8, 5, 2, SequenceMode.BIDIRECTIONAL));
    assertEquals(4, sizeCheck(8, 5, 2, SequenceMode.BIDIRECTIONAL));
    assertEquals(0, sizeCheck(8, 5, 2, SequenceMode.TRANSLATED));
    assertEquals(12, sizeCheck(24, 5, 2, SequenceMode.TRANSLATED));
  }

  private static long sizeCheck(final long length, final int word, final int step, final SequenceMode mode) throws IOException {
    final ReaderParams rp = new MockReaderParams(length, 1L, mode.codeType());
    final ISequenceParams spa = new MockSequenceParams(rp, mode);
    return BuildParams.size(word, step, spa);
  }

  private static final String DNA = ""
          + ">a\n"
          + "acgtacgtacgtacgtacgtacgtacgtacgt\n"
          + ">b\n"
          + "acgtacgtacgtacgtacgtacgtacgtacgt\n"
          + ">c\n"
          + "acgtacgtacgtacgtacgtacgtacgtacgt\n"
          + ">d\n"
          + "acgtacgtacgtacgtacgtacgtacgtacgt\n"
          + ">e\n"
          + "acgtacgtacgtacgtacgtacgtacgtacgt\n";

  public void testSizeB() throws IOException {
    final File dir = FileUtils.createTempDir("test", "buildparamssize");
    try {
      ReaderTestUtils.getDNADir(DNA, dir);
      SequenceParams seq = SequenceParams.builder().directory(dir).readerRestriction(LongRange.NONE).mode(SequenceMode.UNIDIRECTIONAL).create();
      assertEquals(160, seq.reader().lengthBetween(0, 5));
      assertEquals(160, seq.reader().totalLength());
      assertEquals(10, BuildParams.builder().sequences(seq).stepSize(16).windowSize(16).create().calculateSize());
      seq = SequenceParams.builder().directory(dir).readerRestriction(new LongRange(2, 5)).mode(SequenceMode.UNIDIRECTIONAL).useMemReader(false).create();
      assertEquals(3, seq.reader().numberSequences());
      assertEquals(96, seq.reader().lengthBetween(0, 3));
      assertEquals(160, seq.reader().totalLength());
      assertEquals(6, BuildParams.builder().sequences(seq).stepSize(16).windowSize(16).create().calculateSize());
    } finally {
      FileHelper.deleteAll(dir);
    }
  }
}
