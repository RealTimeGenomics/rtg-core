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

package com.rtg.assembler;

import static com.rtg.util.StringUtils.LS;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import com.rtg.mode.DnaUtils;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.reader.SdfId;
import com.rtg.reader.SequencesReader;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.intervals.LongRange;
import com.rtg.util.test.FileHelper;

import junit.framework.TestCase;

/**
 */
public class ReadPairSourceTest extends TestCase {
  static final String LEFT_SEQUENCE = ">a" + LS + "ACGT" + LS + ">b" + LS + "GGGG" + LS;
  static final String RIGHT_SEQUENCE = ">a" + LS + "CCCCAA" + LS + ">b" + LS + "ACAAAC" + LS;

  @Override
  public void setUp() {
    Diagnostic.setLogStream();
  }

  static String fragementToString(List<byte[]> fragments) {
    final StringBuilder sb = new StringBuilder();
    sb.append("[");
    String join = "";
    for (byte[] b : fragments) {
      sb.append(join);
      sb.append(DnaUtils.bytesToSequenceIncCG(b));
      join = ", ";
    }
    sb.append("]");
    return sb.toString();
  }

  static void listEquals(List<byte[]> a, List<byte[]> b) {
    assertEquals("expected <" + fragementToString(a) + "> but was <" + fragementToString(b) + ">", a.size(), b.size());
    for (int i = 0; i < a.size(); ++i) {
      assertTrue("expected <" + fragementToString(a) + "> but was <" + fragementToString(b) + ">", Arrays.equals(a.get(i), b.get(i)));
    }
  }

  public void testReadPairSource() throws IOException {
    final SequencesReader left = ReaderTestUtils.getReaderDnaMemory(LEFT_SEQUENCE);
    final SequencesReader middle = ReaderTestUtils.getReaderDnaMemory(">a" + LS + "AAAA" + LS + ">b" + LS + "ATAT" + LS);
    final SequencesReader right = ReaderTestUtils.getReaderDnaMemory(RIGHT_SEQUENCE);
    final ReadPairSource source = new ReadPairSource(left, middle, right);
    List<byte[]> fragments = source.nextFragments();
    assertNotNull(fragments);
    listEquals(Arrays.asList(DnaUtils.encodeString("ACGT")
        , DnaUtils.encodeString("AAAA")
        , DnaUtils.encodeString("CCCCAA")
    ), fragments);


    fragments = source.nextFragments();
    assertNotNull(fragments);
    listEquals(Arrays.asList(DnaUtils.encodeString("GGGG")
        , DnaUtils.encodeString("ATAT")
        , DnaUtils.encodeString("ACAAAC")
    ), fragments);

    assertNull(source.nextFragments());
    assertNull(source.nextFragments());

    source.reset();
    assertNotNull(source.nextFragments());
    fragments = source.nextFragments();
    assertNotNull(fragments);
    listEquals(Arrays.asList(DnaUtils.encodeString("GGGG")
        , DnaUtils.encodeString("ATAT")
        , DnaUtils.encodeString("ACAAAC")
    ), fragments);

    source.close();
    assertEquals(2, source.numberReads());
  }

  public void testMakeSourceSingle() throws IOException {
    final File tmpDir = FileHelper.createTempDirectory();
    try {
      ReaderTestUtils.getReaderDNA(LEFT_SEQUENCE, tmpDir, new SdfId());
      final ReadPairSource readPairSource = ReadPairSource.makeSource(tmpDir, LongRange.NONE);
      final List<byte[]> fragments = readPairSource.nextFragments();
      assertNotNull(fragments);
      final List<byte[]> expected = new ArrayList<>();
      expected.add(DnaUtils.encodeString("ACGT"));
      listEquals(expected, fragments);
      readPairSource.close();
    } finally {
      FileHelper.deleteAll(tmpDir);
    }
  }

  public void testMakeSourcePaired() throws IOException {
    final File tmpDir = FileHelper.createTempDirectory();
    try {
      ReaderTestUtils.createPairedReaderDNA(LEFT_SEQUENCE, RIGHT_SEQUENCE, tmpDir, new SdfId());
      final ReadPairSource readPairSource = ReadPairSource.makeSource(tmpDir, LongRange.NONE);
      final List<byte[]> fragments = readPairSource.nextFragments();
      assertNotNull(fragments);
      listEquals(Arrays.asList(DnaUtils.encodeString("ACGT")
          , DnaUtils.encodeString("CCCCAA")
      ), fragments);
      readPairSource.close();
    } finally {
      FileHelper.deleteAll(tmpDir);
    }
  }

  public void testExceptions() throws IOException {
    try {
      new ReadPairSource();
      fail();
    } catch (IllegalArgumentException e) {
      // expected
    }
    final SequencesReader left = ReaderTestUtils.getReaderDnaMemory(LEFT_SEQUENCE);
    final SequencesReader empty = ReaderTestUtils.getReaderDnaMemory("");
    try {
      new ReadPairSource(left, empty);
      fail();
    } catch (IllegalArgumentException e) {
      // expected
    }
  }

  public void testSetInsertSizes() throws IOException {
    final SequencesReader left = ReaderTestUtils.getReaderDnaMemory(LEFT_SEQUENCE);
    final SequencesReader right = ReaderTestUtils.getReaderDnaMemory(RIGHT_SEQUENCE);
    final ReadPairSource source = new ReadPairSource(left, right);
    assertEquals(-1, source.maxInsertSize());
    assertEquals(-1, source.minInsertSize());
    source.setMaxInsertSize(10);
    source.setMinInsertSize(4);
    assertEquals(10, source.maxInsertSize());
    assertEquals(4, source.minInsertSize());

  }

  public void testReset() throws IOException {
    // Test that reset will start at first fragment
    final SequencesReader left = ReaderTestUtils.getReaderDnaMemory(LEFT_SEQUENCE);
    final SequencesReader right = ReaderTestUtils.getReaderDnaMemory(RIGHT_SEQUENCE);
    final ReadPairSource source = new ReadPairSource(left, right);
    List<byte[]> fragments = source.nextFragments();
    assertNotNull(fragments);

    final List<byte[]> expected = Arrays.asList(DnaUtils.encodeString("ACGT")
        , DnaUtils.encodeString("CCCCAA")
    );
    listEquals(expected, fragments);
    source.reset();
    fragments = source.nextFragments();
    listEquals(expected, fragments);

  }
}
