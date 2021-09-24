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
