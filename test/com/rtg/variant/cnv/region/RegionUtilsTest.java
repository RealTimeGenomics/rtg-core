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
import com.rtg.util.StringUtils;
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
    assertFalse(r.isInRegion(0));
    assertFalse(r.isInRegion(10));
    assertFalse(r.isInRegion(Integer.MAX_VALUE));
    assertFalse(r.isInRegion(Integer.MIN_VALUE));
    for (int i = 1; i < 10; i++) {
      assertTrue(r.isInRegion(i));
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
    assertFalse(r.isInRegion(0));
    assertFalse(r.isInRegion(3));
    assertFalse(r.isInRegion(6));
    assertFalse(r.isInRegion(10));
    assertFalse(r.isInRegion(Integer.MAX_VALUE));
    assertFalse(r.isInRegion(Integer.MIN_VALUE));
    assertTrue(r.isInRegion(1));
    assertTrue(r.isInRegion(2));
    assertTrue(r.isInRegion(7));
    assertTrue(r.isInRegion(8));
    assertTrue(r.isInRegion(9));
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
    assertFalse(r.isInRegion(0));
    assertFalse(r.isInRegion(3));
    assertFalse(r.isInRegion(4));
    assertFalse(r.isInRegion(5));
    assertFalse(r.isInRegion(6));
    assertFalse(r.isInRegion(10));
    assertFalse(r.isInRegion(Integer.MAX_VALUE));
    assertFalse(r.isInRegion(Integer.MIN_VALUE));
    assertTrue(r.isInRegion(1));
    assertTrue(r.isInRegion(2));
    assertTrue(r.isInRegion(7));
    assertTrue(r.isInRegion(8));
    assertTrue(r.isInRegion(9));

    final Region f = m.get("foo");
    assertNotNull(f);
    assertFalse(f.isInRegion(6));
    assertFalse(f.isInRegion(10));
    assertFalse(f.isInRegion(Integer.MAX_VALUE));
    assertFalse(f.isInRegion(Integer.MIN_VALUE));
    assertTrue(f.isInRegion(7));
    assertTrue(f.isInRegion(8));
    assertTrue(f.isInRegion(9));
}

  private static final String SEQUENCES = ">seq1" + StringUtils.LS
  + "aaaannaaca"
  + "nnnnnnnnnn"
  + "nacacggttc"
  + "ancagttacn"
  + "nncagtcagc"
  + "atnn" + StringUtils.LS

  + ">seq2" + StringUtils.LS
  + "nnn" + StringUtils.LS

  + ">seq3" + StringUtils.LS
  + "nnnactacga"
  + "gcatgact" + StringUtils.LS

  + ">seq4" + StringUtils.LS
  + "acgatcagtc"
  + "agctnnn" + StringUtils.LS

  + ">seq5" + StringUtils.LS
  + "acgtacgtactg" + StringUtils.LS;

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
        for (int i = 0; i <= 3; i++) {
          assertFalse(r.isInRegion(i));
        }
        for (int i = 4; i <= 6; i++) {
          assertTrue(r.isInRegion(i));
        }
        for (int i = 7; i <= 13; i++) {
          assertFalse(r.isInRegion(i));
        }
        assertTrue(r.isInRegion(14));
        for (int i = 15; i <= 20; i++) {
          assertFalse(r.isInRegion(i));
        }
        assertFalse(r.isInRegion(Integer.MAX_VALUE));
        assertFalse(r.isInRegion(Integer.MIN_VALUE));
        r = map.get("seq2");
        assertNotNull(r);
        assertFalse(r.isInRegion(0));
        assertTrue(r.isInRegion(1));
        assertFalse(r.isInRegion(2));
        assertFalse(r.isInRegion(Integer.MAX_VALUE));
        assertFalse(r.isInRegion(Integer.MIN_VALUE));
        r = map.get("seq3");
        assertNotNull(r);
        assertFalse(r.isInRegion(0));
        assertTrue(r.isInRegion(1));
        for (int i = 2; i <= 6; i++) {
          assertFalse(r.isInRegion(i));
        }
        assertFalse(r.isInRegion(Integer.MAX_VALUE));
        assertFalse(r.isInRegion(Integer.MIN_VALUE));
        r = map.get("seq4");
        assertNotNull(r);
        for (int i = 1; i <= 4; i++) {
          assertFalse(r.isInRegion(i));
        }
        assertTrue(r.isInRegion(5));
        assertFalse(r.isInRegion(6));
        assertFalse(r.isInRegion(7));
        assertFalse(r.isInRegion(Integer.MAX_VALUE));
        assertFalse(r.isInRegion(Integer.MIN_VALUE));
        r = map.get("seq5");
        assertNotNull(r);
        assertFalse(r.isInRegion(1));
        assertFalse(r.isInRegion(Integer.MAX_VALUE));
        assertFalse(r.isInRegion(Integer.MIN_VALUE));
        assertEquals(EmptyRegion.EMPTY_REGION, r);

      }
    } finally {
      assertTrue(FileHelper.deleteAll(dir));
    }
  }

  public void testFindGermlineDeletesUnderMean() {
    Region r = RegionUtils.findGermlineDeletesUnderMean(new int[] {10, 10, 10, 0, 0, 10, 10, 10, 0, 0, 0}, 3.0, EmptyRegion.EMPTY_REGION, false);
    for (int i = 1; i <= 3; i++) {
      assertFalse(r.isInRegion(i));
    }
    for (int i = 4; i <= 5; i++) {
      assertTrue(r.isInRegion(i));
    }
    for (int i = 6; i <= 8; i++) {
      assertFalse(r.isInRegion(i));
    }
    for (int i = 9; i <= 11; i++) {
      assertTrue(r.isInRegion(i));
    }
    r = RegionUtils.findGermlineDeletesUnderMean(new int[] {3, 3, 1, 0, 0}, 3.0, EmptyRegion.EMPTY_REGION, false);
    for (int i = 1; i <= 2; i++) {
      assertFalse(r.isInRegion(i));
    }
    for (int i = 3; i <= 5; i++) {
      assertTrue(r.isInRegion(i));
    }
    r = RegionUtils.findGermlineDeletesUnderMean(new int[] {3, 3, 1, 1, 2}, 3.0, EmptyRegion.EMPTY_REGION, false);
    for (int i = 1; i <= 5; i++) {
      assertFalse(r.isInRegion(i));
    }
  }

  public void testFindGermlineDeletesOverMean() {
    Region r = RegionUtils.findGermlineDeletesOverMean(new int[] {9, 9, 9, 9, 9, 100, 100, 5, 5, 5, 5, 5, 5, 5, 5, 100, 100, 100}, 3.0, EmptyRegion.EMPTY_REGION, false);
    for (int i = 1; i <= 5; i++) {
      assertFalse(r.isInRegion(i));
    }
    for (int i = 6; i <= 7; i++) {
      assertTrue(r.isInRegion(i));
    }
    for (int i = 8; i <= 15; i++) {
      assertFalse(r.isInRegion(i));
    }
    for (int i = 16; i <= 18; i++) {
      assertTrue(r.isInRegion(i));
    }
    assertFalse(r.isInRegion(19));
    r = RegionUtils.findGermlineDeletesOverMean(new int[] {9, 9, 9, 9, 9, 100, 100, 5, 5, 5, 5, 5, 5, 5, 5}, 3.0, EmptyRegion.EMPTY_REGION, false);
    for (int i = 1; i <= 4; i++) {
      assertFalse(r.isInRegion(i));
    }
    for (int i = 5; i <= 7; i++) {
      assertTrue(r.isInRegion(i));
    }
    for (int i = 8; i <= 15; i++) {
      assertFalse(r.isInRegion(i));
    }
    r = RegionUtils.findGermlineDeletesOverMean(new int[] {3, 3, 1, 1, 2}, 3.0, EmptyRegion.EMPTY_REGION, false);
    for (int i = 1; i <= 5; i++) {
      assertFalse(r.isInRegion(i));
    }
  }

  public void testNonEmptyIgnoreRegion() {
    final Region r = RegionUtils.findGermlineDeletesUnderMean(new int[] {10, 10, 10, 0, 0, 10, 10, 10, 0, 0, 0}, 3.0, new SimpleCnvRegion(6, 9), false);
    for (int i = 1; i <= 8; i++) {
      assertFalse(r.isInRegion(i));
    }
    for (int i = 9; i <= 11; i++) {
      assertTrue(r.isInRegion(i));
    }
  }
}
