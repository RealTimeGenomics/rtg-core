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
package com.rtg.variant.cnv.region;

import static com.rtg.util.StringUtils.LS;

import java.io.File;
import java.io.IOException;
import java.io.Reader;
import java.io.StringReader;
import java.util.Map;

import com.rtg.reader.ReaderTestUtils;
import com.rtg.reader.SequencesReader;
import com.rtg.reader.SequencesReaderFactory;
import com.rtg.util.TestUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.FileUtils;
import com.rtg.util.test.FileHelper;

import junit.framework.TestCase;

/**
 */
public class RegionUtilsTest extends TestCase {

  @Override
  public void setUp() {
    Diagnostic.setLogStream();
  }

  public void testState() {
    TestUtils.testEnum(RegionUtils.State.class, "[IN, OUT]");
  }

  public void testRegionFromReader0() throws IOException {
    final String inStr = ""
      ;
    final Reader in = new StringReader(inStr);
    final Map<String, Region> m = RegionUtils.regionsFromReader(in);
    assertTrue(m.isEmpty());
  }

  public void testRegionFromReader1() throws IOException {
    final String inStr = ""
      + "name 1 10" + LS
      ;
    final Reader in = new StringReader(inStr);
    final Map<String, Region> m = RegionUtils.regionsFromReader(in);
    assertEquals(1, m.size());
    final Region r = m.get("name");
    assertNotNull(r);
    assertFalse(r.contains(0));
    assertFalse(r.contains(10));
    assertFalse(r.contains(Integer.MAX_VALUE));
    assertFalse(r.contains(Integer.MIN_VALUE));
    for (int i = 1; i < 10; ++i) {
      assertTrue(r.contains(i));
    }
  }

  public void testRegionFromReader2() throws IOException {
    final String inStr = ""
      + "name 1 3" + LS
      + "name 7 10" + LS
      ;
    final Reader in = new StringReader(inStr);
    final Map<String, Region> m = RegionUtils.regionsFromReader(in);
    assertEquals(1, m.size());
    final Region r = m.get("name");
    assertNotNull(r);
    assertFalse(r.contains(0));
    assertFalse(r.contains(3));
    assertFalse(r.contains(6));
    assertFalse(r.contains(10));
    assertFalse(r.contains(Integer.MAX_VALUE));
    assertFalse(r.contains(Integer.MIN_VALUE));
    assertTrue(r.contains(1));
    assertTrue(r.contains(2));
    assertTrue(r.contains(7));
    assertTrue(r.contains(8));
    assertTrue(r.contains(9));
  }

  public void testMultipleNames() throws IOException {
    final String inStr = ""
      + "name 1 3" + LS
      + "foo 7 10" + LS
      + "name 7 10" + LS
      ;
    final Reader in = new StringReader(inStr);
    final Map<String, Region> m = RegionUtils.regionsFromReader(in);
    assertEquals(2, m.size());
    final Region r = m.get("name");
    assertNotNull(r);
    assertFalse(r.contains(0));
    assertFalse(r.contains(3));
    assertFalse(r.contains(4));
    assertFalse(r.contains(5));
    assertFalse(r.contains(6));
    assertFalse(r.contains(10));
    assertFalse(r.contains(Integer.MAX_VALUE));
    assertFalse(r.contains(Integer.MIN_VALUE));
    assertTrue(r.contains(1));
    assertTrue(r.contains(2));
    assertTrue(r.contains(7));
    assertTrue(r.contains(8));
    assertTrue(r.contains(9));

    final Region f = m.get("foo");
    assertNotNull(f);
    assertFalse(f.contains(6));
    assertFalse(f.contains(10));
    assertFalse(f.contains(Integer.MAX_VALUE));
    assertFalse(f.contains(Integer.MIN_VALUE));
    assertTrue(f.contains(7));
    assertTrue(f.contains(8));
    assertTrue(f.contains(9));
}

  private static final String SEQUENCES = ">seq1" + LS
  + "aaaannaaca"
  + "nnnnnnnnnn"
  + "nacacggttc"
  + "ancagttacn"
  + "nncagtcagc"
  + "atnn" + LS

  + ">seq2" + LS
  + "nnn" + LS

  + ">seq3" + LS
  + "nnnactacga"
  + "gcatgact" + LS

  + ">seq4" + LS
  + "acgatcagtc"
  + "agctnnn" + LS

  + ">seq5" + LS
  + "acgtacgtactg" + LS;

  public void testRegionFromSdf() throws IOException {
    final File dir = FileUtils.createTempDir("nblockdetect", "test");
    try {
      final File sdf = new File(dir, "sdf");
      ReaderTestUtils.getReaderDNA(SEQUENCES, sdf, null).close();
      try (SequencesReader reader = SequencesReaderFactory.createDefaultSequencesReader(sdf)) {
        final Map<String, Region> map = RegionUtils.regionsFromSDF(reader, 3);
        assertEquals(5, map.size());
        Region r = map.get("seq1");
        assertNotNull(r);
        for (int i = 0; i <= 3; ++i) {
          assertFalse(r.contains(i));
        }
        for (int i = 4; i <= 6; ++i) {
          assertTrue(r.contains(i));
        }
        for (int i = 7; i <= 13; ++i) {
          assertFalse(r.contains(i));
        }
        assertTrue(r.contains(14));
        for (int i = 15; i <= 20; ++i) {
          assertFalse(r.contains(i));
        }
        assertFalse(r.contains(Integer.MAX_VALUE));
        assertFalse(r.contains(Integer.MIN_VALUE));
        r = map.get("seq2");
        assertNotNull(r);
        assertFalse(r.contains(0));
        assertTrue(r.contains(1));
        assertFalse(r.contains(2));
        assertFalse(r.contains(Integer.MAX_VALUE));
        assertFalse(r.contains(Integer.MIN_VALUE));
        r = map.get("seq3");
        assertNotNull(r);
        assertFalse(r.contains(0));
        assertTrue(r.contains(1));
        for (int i = 2; i <= 6; ++i) {
          assertFalse(r.contains(i));
        }
        assertFalse(r.contains(Integer.MAX_VALUE));
        assertFalse(r.contains(Integer.MIN_VALUE));
        r = map.get("seq4");
        assertNotNull(r);
        for (int i = 1; i <= 4; ++i) {
          assertFalse(r.contains(i));
        }
        assertTrue(r.contains(5));
        assertFalse(r.contains(6));
        assertFalse(r.contains(7));
        assertFalse(r.contains(Integer.MAX_VALUE));
        assertFalse(r.contains(Integer.MIN_VALUE));
        r = map.get("seq5");
        assertNotNull(r);
        assertFalse(r.contains(1));
        assertFalse(r.contains(Integer.MAX_VALUE));
        assertFalse(r.contains(Integer.MIN_VALUE));
        assertEquals(EmptyRegion.EMPTY_REGION, r);

      }
    } finally {
      assertTrue(FileHelper.deleteAll(dir));
    }
  }

  public void testFindGermlineDeletesUnderMean() {
    Region r = RegionUtils.findGermlineDeletesUnderMean(new int[] {10, 10, 10, 0, 0, 10, 10, 10, 0, 0, 0}, 3.0, EmptyRegion.EMPTY_REGION, false);
    for (int i = 1; i <= 3; ++i) {
      assertFalse(r.contains(i));
    }
    for (int i = 4; i <= 5; ++i) {
      assertTrue(r.contains(i));
    }
    for (int i = 6; i <= 8; ++i) {
      assertFalse(r.contains(i));
    }
    for (int i = 9; i <= 11; ++i) {
      assertTrue(r.contains(i));
    }
    r = RegionUtils.findGermlineDeletesUnderMean(new int[] {3, 3, 1, 0, 0}, 3.0, EmptyRegion.EMPTY_REGION, false);
    for (int i = 1; i <= 2; ++i) {
      assertFalse(r.contains(i));
    }
    for (int i = 3; i <= 5; ++i) {
      assertTrue(r.contains(i));
    }
    r = RegionUtils.findGermlineDeletesUnderMean(new int[] {3, 3, 1, 1, 2}, 3.0, EmptyRegion.EMPTY_REGION, false);
    for (int i = 1; i <= 5; ++i) {
      assertFalse(r.contains(i));
    }
  }

  public void testFindGermlineDeletesOverMean() {
    Region r = RegionUtils.findGermlineDeletesOverMean(new int[] {9, 9, 9, 9, 9, 100, 100, 5, 5, 5, 5, 5, 5, 5, 5, 100, 100, 100}, 3.0, EmptyRegion.EMPTY_REGION, false);
    for (int i = 1; i <= 5; ++i) {
      assertFalse(r.contains(i));
    }
    for (int i = 6; i <= 7; ++i) {
      assertTrue(r.contains(i));
    }
    for (int i = 8; i <= 15; ++i) {
      assertFalse(r.contains(i));
    }
    for (int i = 16; i <= 18; ++i) {
      assertTrue(r.contains(i));
    }
    assertFalse(r.contains(19));
    r = RegionUtils.findGermlineDeletesOverMean(new int[] {9, 9, 9, 9, 9, 100, 100, 5, 5, 5, 5, 5, 5, 5, 5}, 3.0, EmptyRegion.EMPTY_REGION, false);
    for (int i = 1; i <= 4; ++i) {
      assertFalse(r.contains(i));
    }
    for (int i = 5; i <= 7; ++i) {
      assertTrue(r.contains(i));
    }
    for (int i = 8; i <= 15; ++i) {
      assertFalse(r.contains(i));
    }
    r = RegionUtils.findGermlineDeletesOverMean(new int[] {3, 3, 1, 1, 2}, 3.0, EmptyRegion.EMPTY_REGION, false);
    for (int i = 1; i <= 5; ++i) {
      assertFalse(r.contains(i));
    }
  }

  public void testNonEmptyIgnoreRegion() {
    final Region r = RegionUtils.findGermlineDeletesUnderMean(new int[] {10, 10, 10, 0, 0, 10, 10, 10, 0, 0, 0}, 3.0, new SimpleCnvRegion(6, 9), false);
    for (int i = 1; i <= 8; ++i) {
      assertFalse(r.contains(i));
    }
    for (int i = 9; i <= 11; ++i) {
      assertTrue(r.contains(i));
    }
  }
}
