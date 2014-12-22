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
package com.rtg.variant.bayes.multisample;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Deque;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import com.rtg.mode.DNA;
import com.rtg.variant.Variant;
import com.rtg.variant.bayes.multisample.population.AlleleCounts;
import com.rtg.variant.bayes.multisample.population.SiteSpecificPriors;

import junit.framework.TestCase;

/**
 * Test class
 */
public class ComplexitiesTest extends TestCase {
  /**
   * Generate a dna byte array of a certain length consisting of ACGT repeatedly
   * @param length length of template to generate
   * @return a byte array of the desired length with 0-4 in each position
   */
  public static byte[] template(int length) {
    final byte[] res = new byte[length];
    for (int i = 0; i < length; i++) {
      res[i] = (byte) (i % 5);
    }
    return res;
  }

  public void testSimple() {
    final ArrayList<Variant> chunk = new ArrayList<>();
    chunk.add(TestUtils.createVariant(4));
    chunk.add(TestUtils.createVariant(5));
    chunk.add(TestUtils.createVariant(6));
    chunk.add(TestUtils.createVariant(10));
    chunk.add(TestUtils.createVariant(15));

    final Complexities regions = new Complexities(chunk, "foo", 0, 50, 3, 15, template(30), true, null);
    regions.globalIntegrity();
    assertEquals(2, regions.size());
    assertNotNull(regions.danglingStart());
    assertNotNull(regions.danglingEnd());
    try {
      regions.iterator();
      fail();
    } catch (final IllegalStateException e) {
      assertEquals("Dangling boundaries haven't been fixed (left: false right: false)", e.getMessage());
    }
    final Iterator<ComplexRegion> it = regions.iteratorInternal();
    assertTrue(it.hasNext());
    final ComplexRegion com = it.next();
    assertEquals(4, com.getStart());
    assertEquals(7, com.getEnd());
    assertEquals(ComplexRegion.RegionType.COMPLEX, com.type());
    assertTrue(it.hasNext());

    assertEquals("*foo[0..50)*#5[foo[4..7)CX, foo[10..16)CX]", regions.toString());
  }

  abstract static class MockSsp implements SiteSpecificPriors {

    @Override
    public HaploidDiploidHypotheses<?> getSnpHypotheses(String templateName, int zeroPos) {
      return null;
    }
  }

  public void testSiteSpecificPriors() {
    final ArrayList<Variant> chunk = new ArrayList<>();
    chunk.add(TestUtils.createVariant(4));
    chunk.add(TestUtils.createVariant(5));
    chunk.add(TestUtils.createVariant(6));
    chunk.add(TestUtils.createVariant(10));
    chunk.add(TestUtils.createVariant(15));
    final Map<String, Integer> empty = new HashMap<>();
    Complexities regions = new Complexities(chunk, "foo", 0, 50, 3, 15, template(30), false, new MockSsp() {
      @Override
      public List<AlleleCounts> getCounts(final String templateName, int start, int end) {
        return Arrays.asList(new AlleleCounts(2, empty, "A"));
      }
    });
    regions.globalIntegrity();
    assertEquals(1, regions.size());

    Iterator<ComplexRegion> it = regions.iteratorInternal();
    assertTrue(it.hasNext());
    ComplexRegion com = it.next();
    assertEquals(2, com.getStart());
    assertEquals(7, com.getEnd());
    assertEquals(ComplexRegion.RegionType.COMPLEX, com.type());
    assertFalse(it.hasNext());

    regions = new Complexities(chunk, "foo", 0, 50, 3, 15, template(30), false, new MockSsp() {
      @Override
      public List<AlleleCounts> getCounts(final String templateName, int start, int end) {
        return Arrays.asList(new AlleleCounts(7, empty, "C"));
      }
    });
    regions.globalIntegrity();
    assertEquals(1, regions.size());

    it = regions.iteratorInternal();
    assertTrue(it.hasNext());
    com = it.next();
    assertEquals(4, com.getStart());
    assertEquals(11, com.getEnd());
    assertEquals(ComplexRegion.RegionType.COMPLEX, com.type());
    assertFalse(it.hasNext());

    regions = new Complexities(chunk, "foo", 0, 50, 3, 15, template(30), false, new MockSsp() {
      @Override
      public List<AlleleCounts> getCounts(final String templateName, int start, int end) {
        return Arrays.asList(new AlleleCounts(12, empty, "G"));
      }
    });
    regions.globalIntegrity();
    assertEquals(2, regions.size());

    it = regions.iteratorInternal();
    assertTrue(it.hasNext());
    com = it.next();
    assertEquals(4, com.getStart());
    assertEquals(7, com.getEnd());
    assertEquals(ComplexRegion.RegionType.COMPLEX, com.type());

    assertTrue(it.hasNext());
    com = it.next();
    assertEquals(10, com.getStart());
    assertEquals(16, com.getEnd());
    assertEquals(ComplexRegion.RegionType.COMPLEX, com.type());
    assertFalse(it.hasNext());

    regions = new Complexities(chunk, "foo", 0, 50, 3, 15, template(30), false, new MockSsp() {
      @Override
      public List<AlleleCounts> getCounts(final String templateName, int start, int end) {
        return Arrays.asList(new AlleleCounts(7, empty, "T"), new AlleleCounts(12, empty, "A"));
      }
    });
    regions.globalIntegrity();
    assertEquals(1, regions.size());

    it = regions.iteratorInternal();
    assertTrue(it.hasNext());
    com = it.next();
    assertEquals(4, com.getStart());
    assertEquals(16, com.getEnd());
    assertEquals(ComplexRegion.RegionType.COMPLEX, com.type());

    regions = new Complexities(chunk, "foo", 0, 50, 3, 15, template(30), false, new MockSsp() {
      @Override
      public List<AlleleCounts> getCounts(final String templateName, int start, int end) {
        return Arrays.asList(new AlleleCounts(7, empty, "ACGTACGT"));
      }
    });
    regions.globalIntegrity();
    assertEquals(1, regions.size());

    it = regions.iteratorInternal();
    assertTrue(it.hasNext());
    com = it.next();
    assertEquals(4, com.getStart());
    assertEquals(16, com.getEnd());
    assertEquals(ComplexRegion.RegionType.COMPLEX, com.type());
  }

  //a region that is not interesting (shouldn't set anything).
  public void testBoringRegion() {
    final ArrayList<Variant> chunk = new ArrayList<>();
    chunk.add(TestUtils.createVariant(3, 4, false, 0, false));
    chunk.add(TestUtils.createVariant(4));
    chunk.add(TestUtils.createVariant(5));
    chunk.add(TestUtils.createVariant(6));
    chunk.add(TestUtils.createVariant(10));
    chunk.add(TestUtils.createVariant(15));

    final Complexities regions = new Complexities(chunk, "foo", 0, 50, 3, 15, template(30), false, null);
    regions.globalIntegrity();
    assertEquals(1, regions.size());
    assertEquals(4, regions.danglingStart().getStart());
    assertEquals(15, regions.danglingEnd().getStart());
    try {
      regions.iterator();
      fail();
    } catch (final IllegalStateException e) {
      assertEquals("Dangling boundaries haven't been fixed (left: false right: false)", e.getMessage());
    }
    final Iterator<ComplexRegion> it = regions.iteratorInternal();
    assertTrue(it.hasNext());
    final ComplexRegion com = it.next();
    assertEquals(4, com.getStart());
    assertEquals(7, com.getEnd());
    assertEquals(ComplexRegion.RegionType.COMPLEX, com.type());
    assertFalse(it.hasNext());
  }

  public void testForceComplex() {
    final ArrayList<Variant> chunk = new ArrayList<>();
    chunk.add(TestUtils.createVariant(4));
    chunk.add(TestUtils.createVariant(5));
    chunk.add(TestUtils.createVariant(6));
    chunk.add(TestUtils.createVariant(10, true));
    chunk.add(TestUtils.createVariant(15));

    final Complexities regions = new Complexities(chunk, "foo", 0, 50, 3, 15, template(30), false, null);
    regions.globalIntegrity();
    assertEquals(2, regions.size());
    final Iterator<ComplexRegion> it = regions.iteratorInternal();
    assertTrue(it.hasNext());
    final ComplexRegion com = it.next();
    assertEquals(4, com.getStart());
    assertEquals(7, com.getEnd());
    assertEquals(ComplexRegion.RegionType.COMPLEX, com.type());
    assertTrue(it.hasNext());
    final ComplexRegion com2 = it.next();
    assertEquals(10, com2.getStart());
    assertEquals(11, com2.getEnd());
    assertEquals(ComplexRegion.RegionType.COMPLEX, com2.type());
    assertFalse(it.hasNext());
  }

  public void testDangling() {
    final ArrayList<Variant> chunk = new ArrayList<>();
    chunk.add(TestUtils.createVariant(4));
    chunk.add(TestUtils.createVariant(5));
    chunk.add(TestUtils.createVariant(6));
    chunk.add(TestUtils.createVariant(10, true));

    final Complexities regions = new Complexities(chunk, "foo", 1, 14, 3, 10, template(30), true, null);
    regions.globalIntegrity();
    assertEquals(2, regions.size());
    assertNotNull(regions.danglingStart());
    assertEquals(4, regions.danglingStart().getStart());
    assertEquals(10, regions.danglingEnd().getStart());
    final Iterator<ComplexRegion> it = regions.iteratorInternal();
    assertTrue(it.hasNext());
    final ComplexRegion com = it.next();
    assertEquals(4, com.getStart());
    assertEquals(7, com.getEnd());
    assertEquals(ComplexRegion.RegionType.COMPLEX, com.type());
    assertTrue(it.hasNext());
    final ComplexRegion com2 = it.next();
    assertEquals(10, com2.getStart());
    assertEquals(11, com2.getEnd());
    assertEquals(ComplexRegion.RegionType.COMPLEX, com2.type());
    assertFalse(it.hasNext());
  }

  public void testFirstPass() {
    final Complexities complex = new Complexities(Collections.<Variant>emptyList(), "foo", 0, 50, 3, 15, template(30), false, null);
    complex.globalIntegrity();
    assertEquals(0, complex.size());
    final ArrayList<Variant> chunk = new ArrayList<>();
    chunk.add(TestUtils.createVariant(0));
    chunk.add(TestUtils.createVariant(4));
    chunk.add(TestUtils.createVariant(5));
    chunk.add(TestUtils.createVariant(6));
    chunk.add(TestUtils.createVariant(10, true));
    chunk.add(TestUtils.createVariant(15));
    final Deque<ComplexRegion> regions = complex.getComplexRegions(chunk);
    final Iterator<ComplexRegion> it = regions.iterator();
    assertTrue(it.hasNext());
    final ComplexRegion com0 = it.next();
    assertEquals(0, com0.getStart());
    assertEquals(1, com0.getEnd());
    assertEquals(ComplexRegion.RegionType.INTERESTING, com0.type());
    assertTrue(it.hasNext());
    final ComplexRegion com1 = it.next();
    assertEquals(4, com1.getStart());
    assertEquals(7, com1.getEnd());
    assertEquals(ComplexRegion.RegionType.COMPLEX, com1.type());
    assertTrue(it.hasNext());
    final ComplexRegion com2 = it.next();
    assertEquals(10, com2.getStart());
    assertEquals(11, com2.getEnd());
    assertEquals(ComplexRegion.RegionType.COMPLEX, com2.type());
    assertTrue(it.hasNext());
    final ComplexRegion com3 = it.next();
    assertEquals(15, com3.getStart());
    assertEquals(16, com3.getEnd());
    assertEquals(ComplexRegion.RegionType.INTERESTING, com3.type());
    assertFalse(it.hasNext());
  }

  public void testIndelExtension() {
    final Complexities complex = new Complexities(Collections.<Variant>emptyList(), "foo", 0, 50, 3, 15, template(30), false, null);
    complex.globalIntegrity();
    assertEquals(0, complex.size());
    final ArrayList<Variant> chunk = new ArrayList<>();
    chunk.add(TestUtils.createVariant(0));
    chunk.add(TestUtils.createVariant(4));
    chunk.add(TestUtils.createVariant(10, 11, true, 3)); // interesting sep + indel length lets it join to previous
    chunk.add(TestUtils.createVariant(15));
    final Deque<ComplexRegion> regions = complex.getComplexRegions(chunk);
    final Iterator<ComplexRegion> it = regions.iterator();
    assertTrue(it.hasNext());
    final ComplexRegion com0 = it.next();
    assertEquals(0, com0.getStart());
    assertEquals(1, com0.getEnd());
    assertEquals(ComplexRegion.RegionType.INTERESTING, com0.type());
    assertTrue(it.hasNext());
    final ComplexRegion com1 = it.next();
    assertEquals(4, com1.getStart());
    assertEquals(11, com1.getEnd());
    assertEquals(ComplexRegion.RegionType.COMPLEX, com1.type());
    assertTrue(it.hasNext());
    final ComplexRegion com3 = it.next();
    assertEquals(15, com3.getStart());
    assertEquals(16, com3.getEnd());
    assertEquals(ComplexRegion.RegionType.INTERESTING, com3.type());
    assertFalse(it.hasNext());
  }

  public void testFirstPassHyper() {
    final Complexities complex = new Complexities(Collections.<Variant>emptyList(), "foo", 0, 50, 3, 5, template(30), false, null);
    assertEquals(0, complex.size());
    final ArrayList<Variant> chunk = new ArrayList<>();
    chunk.add(TestUtils.createVariant(0));
    chunk.add(TestUtils.createVariant(4));
    chunk.add(TestUtils.createVariant(5));
    chunk.add(TestUtils.createVariant(6));
    chunk.add(TestUtils.createVariant(8));
    chunk.add(TestUtils.createVariant(10, true));
    chunk.add(TestUtils.createVariant(15));
    final Deque<ComplexRegion> regions = complex.getComplexRegions(chunk);
    final Iterator<ComplexRegion> it = regions.iterator();
    assertTrue(it.hasNext());
    final ComplexRegion com0 = it.next();
    assertEquals(0, com0.getStart());
    assertEquals(1, com0.getEnd());
    assertEquals(ComplexRegion.RegionType.INTERESTING, com0.type());
    assertTrue(it.hasNext());
    final ComplexRegion com1 = it.next();
    assertEquals(4, com1.getStart());
    assertEquals(11, com1.getEnd());
    assertEquals(ComplexRegion.RegionType.HYPER, com1.type());
    assertTrue(it.hasNext());
    final ComplexRegion com3 = it.next();
    assertEquals(15, com3.getStart());
    assertEquals(16, com3.getEnd());
    assertEquals(ComplexRegion.RegionType.INTERESTING, com3.type());
    assertFalse(it.hasNext());
  }

  public void testStartDangling() {
    final Complexities complex = new Complexities(Collections.<Variant>emptyList(), "foo", 0, 50, 3, 5, template(50), true, null);
    complex.globalIntegrity();
    for (int i = 0; i < 2; i++) {
      final LinkedList<ComplexRegion> regions = new LinkedList<>();
      final ComplexRegion complexRegion = new ComplexRegion("foo", i, i + 1, ComplexRegion.RegionType.INTERESTING);
      regions.add(complexRegion);
      assertTrue("" + i + "expected <" + complexRegion + "> but was <" + complex.computeStartDangle(regions) + ">" , complexRegion == complex.computeStartDangle(regions));
    }
    final LinkedList<ComplexRegion> regions = new LinkedList<>();
    // Cannot dangle off the start of the template so position 4 is the first nong dangling position
    final ComplexRegion complexRegion = new ComplexRegion("foo", 4, 5, ComplexRegion.RegionType.INTERESTING);
    regions.add(complexRegion);
    final ComplexRegion startDangle = complex.computeStartDangle(regions);
    assertNotNull(startDangle);
    assertEquals(4, startDangle.getStart());
    assertNull(complex.computeStartDangle(new LinkedList<ComplexRegion>()));
  }

  public void testEndDangling() {
    final Complexities complex = new Complexities(Collections.<Variant>emptyList(), "foo", 0, 50, 3, 5, template(30), true, null);
    complex.globalIntegrity();
    for (int i = 0; i < 2; i++) {
      final LinkedList<ComplexRegion> regions = new LinkedList<>();
      final ComplexRegion complexRegion = new ComplexRegion("foo", 50 - i - 1, 50 - i, ComplexRegion.RegionType.INTERESTING);
      regions.add(complexRegion);
      assertTrue("" + i, complexRegion == complex.computeEndDangle(regions));
    }
    final LinkedList<ComplexRegion> regions = new LinkedList<>();
    final ComplexRegion complexRegion = new ComplexRegion("foo", 46, 47, ComplexRegion.RegionType.INTERESTING);
    regions.add(complexRegion);
    final ComplexRegion endDangle = complex.computeEndDangle(regions);
    assertNotNull(endDangle);
    assertEquals(46, endDangle.getStart());
    assertNull(complex.computeEndDangle(new LinkedList<ComplexRegion>()));
  }

  public void testMerge() {
    final ComplexRegion a = new ComplexRegion("foo", 4, 7, ComplexRegion.RegionType.INTERESTING);
    final ComplexRegion b = new ComplexRegion("foo", 9, 15, ComplexRegion.RegionType.INTERESTING);
    final ComplexRegion m = Complexities.mergeRegions(a, b, 15);
    assertEquals(4, m.getStart());
    assertEquals(15, m.getEnd());
    assertEquals(ComplexRegion.RegionType.COMPLEX, m.type());
    final ComplexRegion m2 = Complexities.mergeRegions(a, b, 10);
    assertEquals(4, m2.getStart());
    assertEquals(15, m2.getEnd());
    assertEquals(ComplexRegion.RegionType.HYPER, m2.type());
  }

  //no boundary case
  public void testFixDanglingTrivial() {
    final ArrayList<Variant> chunkA = new ArrayList<>();
    chunkA.add(TestUtils.createVariant(7, true));
    final Complexities regionsA = new Complexities(chunkA, "foo", 0, 20, 3, 15, template(30), false, null);
    regionsA.globalIntegrity();
    final ArrayList<Variant> chunkB = new ArrayList<>();
    chunkB.add(TestUtils.createVariant(24, true));
    final Complexities regionsB = new Complexities(chunkB, "foo", 20, 40, 3, 15, template(30), false, null);
    regionsB.globalIntegrity();
    Complexities.fixDangling(regionsA, regionsB);
    regionsA.globalIntegrity();
    regionsB.globalIntegrity();
    assertNotNull(regionsA.danglingStart());
    assertEquals(7, regionsA.danglingStart().getStart());
    assertNull(regionsB.danglingStart());
    assertNull(regionsA.danglingEnd());
    assertNotNull(regionsB.danglingEnd());
    assertEquals(24, regionsB.danglingEnd().getStart());

    assertEquals("*foo[0..20)#1[foo[7..8)CX]", regionsA.toString());
    assertEquals("foo[20..40)*#1[foo[24..25)CX]", regionsB.toString());

  }

  //simple adjacent complex calls on boundary
  public void testFixDangling() {
    final ArrayList<Variant> chunkA = new ArrayList<>();
    chunkA.add(TestUtils.createVariant(19, true));
    final Complexities regionsA = new Complexities(chunkA, "foo", 0, 20, 3, 15, template(30), false, null);
    regionsA.globalIntegrity();
    final ArrayList<Variant> chunkB = new ArrayList<>();
    chunkB.add(TestUtils.createVariant(21, true));
    final Complexities regionsB = new Complexities(chunkB, "foo", 20, 40, 3, 15, template(30), false, null);
    regionsB.globalIntegrity();
    assertNotNull(regionsA.danglingEnd());
    assertEquals(19, regionsA.danglingEnd().getStart());
    assertEquals(20, regionsA.danglingEnd().getEnd());
    assertEquals(ComplexRegion.RegionType.COMPLEX, regionsA.danglingEnd().type());
    assertNotNull(regionsB.danglingStart());
    assertEquals(21, regionsB.danglingStart().getStart());
    assertEquals(22, regionsB.danglingStart().getEnd());
    assertEquals(ComplexRegion.RegionType.COMPLEX, regionsB.danglingStart().type());
    Complexities.fixDangling(null, regionsA);
    Complexities.fixDangling(regionsA, regionsB);
    Complexities.fixDangling(regionsB, null);
    regionsA.globalIntegrity();
    regionsB.globalIntegrity();
    assertNull(regionsA.danglingStart());
    assertNull(regionsB.danglingStart());
    assertNull(regionsA.danglingEnd());
    assertNull(regionsB.danglingEnd());
    assertEquals(1, regionsA.size());
    assertEquals(1, regionsB.size());
    final Iterator<ComplexRegion> itA = regionsA.iterator();
    assertTrue(itA.hasNext());
    final ComplexRegion mergedA = itA.next();
    assertEquals(19, mergedA.getStart());
    assertEquals(22, mergedA.getEnd());
    assertEquals(ComplexRegion.RegionType.COMPLEX, mergedA.type());
    assertFalse(itA.hasNext());

    final Iterator<ComplexRegion> itB = regionsB.iterator();
    assertTrue(itB.hasNext());
    final ComplexRegion mergedB = itB.next();
    assertEquals(19, mergedB.getStart());
    assertEquals(22, mergedB.getEnd());
    assertEquals(ComplexRegion.RegionType.COMPLEX, mergedB.type());
    assertFalse(itB.hasNext());
  }

  //simple adjacent interesting calls on boundary - close enough to merge
  public void testFixDanglingInterestingA() {
    final ArrayList<Variant> chunkA = new ArrayList<>();
    chunkA.add(TestUtils.createVariant(19));
    final Complexities regionsA = new Complexities(chunkA, "foo", 0, 20, 3, 15, template(30), true, null);
    regionsA.globalIntegrity();
    final ArrayList<Variant> chunkB = new ArrayList<>();
    chunkB.add(TestUtils.createVariant(21));
    final Complexities regionsB = new Complexities(chunkB, "foo", 20, 40, 3, 15, template(30), true, null);
    regionsB.globalIntegrity();
    assertNotNull(regionsA.danglingEnd());
    assertEquals(19, regionsA.danglingEnd().getStart());
    assertEquals(20, regionsA.danglingEnd().getEnd());
    assertEquals(ComplexRegion.RegionType.INTERESTING, regionsA.danglingEnd().type());
    assertNotNull(regionsB.danglingStart());
    assertEquals(21, regionsB.danglingStart().getStart());
    assertEquals(22, regionsB.danglingStart().getEnd());
    assertEquals(ComplexRegion.RegionType.INTERESTING, regionsB.danglingStart().type());
    assertFalse(regionsA.isFixed());
    Complexities.fixDangling(null, regionsA);
    Complexities.fixDangling(regionsA, regionsB);
    Complexities.fixDangling(regionsB, null);
    assertTrue(regionsA.isFixed());
    regionsA.globalIntegrity();
    regionsB.globalIntegrity();
    assertNull(regionsA.danglingStart());
    assertNull(regionsB.danglingStart());
    assertNull(regionsA.danglingEnd());
    assertNull(regionsB.danglingEnd());
    assertEquals(1, regionsA.size());
    assertEquals(1, regionsB.size());
    final Iterator<ComplexRegion> it = regionsA.iterator();
    assertTrue(it.hasNext());
    final ComplexRegion merged = it.next();
    assertEquals(19, merged.getStart());
    assertEquals(22, merged.getEnd());
    assertEquals(ComplexRegion.RegionType.COMPLEX, merged.type());
    assertFalse(it.hasNext());

    final Iterator<ComplexRegion> itB = regionsB.iterator();
    assertTrue(itB.hasNext());
    final ComplexRegion mergedB = itB.next();
    assertEquals(19, mergedB.getStart());
    assertEquals(22, mergedB.getEnd());
    assertEquals(ComplexRegion.RegionType.COMPLEX, mergedB.type());
    assertFalse(itB.hasNext());

    //check that single chunk does the same thing
    final ArrayList<Variant> chunkC = new ArrayList<>();
    chunkC.add(TestUtils.createVariant(19));
    chunkB.add(TestUtils.createVariant(21));
    final Complexities regionsC = new Complexities(chunkC, "foo", 20, 40, 3, 15, template(30), true, null);
    regionsC.globalIntegrity();
    final Iterator<ComplexRegion> itC = regionsB.iterator();
    assertTrue(itC.hasNext());
    final ComplexRegion mergedC = itC.next();
    assertEquals(19, mergedC.getStart());
    assertEquals(22, mergedC.getEnd());
    assertEquals(ComplexRegion.RegionType.COMPLEX, mergedC.type());
    assertFalse(itC.hasNext());
  }

  //simple adjacent interesting calls on boundary - one further apart still merge
  //Bug 1399
  public void testFixDanglingInterestingB() {
    final ArrayList<Variant> chunkA = new ArrayList<>();
    chunkA.add(TestUtils.createVariant(19));
    final Complexities regionsA = new Complexities(chunkA, "foo", 0, 20, 3, 15, template(30), true, null);
    regionsA.globalIntegrity();
    final ArrayList<Variant> chunkB = new ArrayList<>();
    chunkB.add(TestUtils.createVariant(22));
    final Complexities regionsB = new Complexities(chunkB, "foo", 20, 40, 3, 15, template(30), true, null);
    regionsB.globalIntegrity();
    assertNotNull(regionsA.danglingEnd());
    assertEquals(19, regionsA.danglingEnd().getStart());
    assertEquals(20, regionsA.danglingEnd().getEnd());
    assertEquals(ComplexRegion.RegionType.INTERESTING, regionsA.danglingEnd().type());
    assertNotNull(regionsB.danglingStart());
    assertEquals(22, regionsB.danglingStart().getStart());
    assertEquals(23, regionsB.danglingStart().getEnd());
    assertEquals(ComplexRegion.RegionType.INTERESTING, regionsB.danglingStart().type());
    Complexities.fixDangling(null, regionsA);
    Complexities.fixDangling(regionsA, regionsB);
    Complexities.fixDangling(regionsB, null);
    regionsA.globalIntegrity();
    regionsB.globalIntegrity();
    assertNull(regionsA.danglingStart());
    assertNull(regionsB.danglingStart());
    assertNull(regionsA.danglingEnd());
    assertNull(regionsB.danglingEnd());
    assertEquals(1, regionsA.size());
    assertEquals(1, regionsB.size());
    final Iterator<ComplexRegion> it = regionsA.iterator();
    assertTrue(it.hasNext());
    final ComplexRegion merged = it.next();
    assertEquals(19, merged.getStart());
    assertEquals(23, merged.getEnd());
    assertEquals(ComplexRegion.RegionType.COMPLEX, merged.type());
    assertFalse(it.hasNext());

    final Iterator<ComplexRegion> itB = regionsB.iterator();
    assertTrue(itB.hasNext());
    final ComplexRegion mergedB = itB.next();
    assertEquals(19, mergedB.getStart());
    assertEquals(23, mergedB.getEnd());
    assertEquals(ComplexRegion.RegionType.COMPLEX, mergedB.type());
    assertFalse(itB.hasNext());

    //check that single chunk does the same thing
    final ArrayList<Variant> chunkC = new ArrayList<>();
    chunkC.add(TestUtils.createVariant(19));
    chunkB.add(TestUtils.createVariant(21));
    final Complexities regionsC = new Complexities(chunkC, "foo", 20, 40, 3, 15, template(30), true, null);
    regionsC.globalIntegrity();
    final Iterator<ComplexRegion> itC = regionsB.iterator();
    assertTrue(itC.hasNext());
    final ComplexRegion mergedC = itC.next();
    assertEquals(19, mergedC.getStart());
    assertEquals(23, mergedC.getEnd());
    assertEquals(ComplexRegion.RegionType.COMPLEX, mergedC.type());
    assertFalse(itC.hasNext());
  }

  //simple adjacent interesting calls on boundary - just far apart enough not to merge (still generate dangling though)
  //Bug 1399
  public void testFixDanglingInterestingC() {
    final ArrayList<Variant> chunkA = new ArrayList<>();
    chunkA.add(TestUtils.createVariant(18));
    final Complexities regionsA = new Complexities(chunkA, "foo", 0, 20, 3, 15, template(30), true, null);
    regionsA.globalIntegrity();
    final ArrayList<Variant> chunkB = new ArrayList<>();
    chunkB.add(TestUtils.createVariant(22));
    final Complexities regionsB = new Complexities(chunkB, "foo", 20, 40, 3, 15, template(30), true, null);
    regionsB.globalIntegrity();
    assertNotNull(regionsA.danglingEnd());
    assertEquals(18, regionsA.danglingEnd().getStart());
    assertEquals(19, regionsA.danglingEnd().getEnd());
    assertEquals(ComplexRegion.RegionType.INTERESTING, regionsA.danglingEnd().type());
    assertNotNull(regionsB.danglingStart());
    assertEquals(22, regionsB.danglingStart().getStart());
    assertEquals(23, regionsB.danglingStart().getEnd());
    assertEquals(ComplexRegion.RegionType.INTERESTING, regionsB.danglingStart().type());
    Complexities.fixDangling(null, regionsA);
    Complexities.fixDangling(regionsA, regionsB);
    Complexities.fixDangling(regionsB, null);
    regionsA.globalIntegrity();
    regionsB.globalIntegrity();
    assertNull(regionsA.danglingStart());
    assertNull(regionsB.danglingStart());
    assertNull(regionsA.danglingEnd());
    assertNull(regionsB.danglingEnd());
    assertEquals(0, regionsA.size());
    assertEquals(0, regionsB.size());
    final Iterator<ComplexRegion> it = regionsA.iterator();
    assertFalse(it.hasNext());

    final Iterator<ComplexRegion> itB = regionsB.iterator();
    assertFalse(itB.hasNext());

    //check that single chunk does the same thing
    final ArrayList<Variant> chunkC = new ArrayList<>();
    chunkC.add(TestUtils.createVariant(18));
    chunkB.add(TestUtils.createVariant(21));
    final Complexities regionsC = new Complexities(chunkC, "foo", 20, 40, 3, 15, template(30), true, null);
    regionsC.globalIntegrity();
    final Iterator<ComplexRegion> itC = regionsB.iterator();
    assertFalse(itC.hasNext());
  }

  //simple adjacent interesting calls on boundary - the second is right on the start of the second chunk (more boundary cases)
  public void testFixDanglingInterestingx() {
    final ArrayList<Variant> chunkA = new ArrayList<>();
    chunkA.add(TestUtils.createVariant(18));
    final Complexities regionsA = new Complexities(chunkA, "foo", 0, 20, 3, 15, template(30), true, null);
    regionsA.globalIntegrity();
    final ArrayList<Variant> chunkB = new ArrayList<>();
    chunkB.add(TestUtils.createVariant(20));
    final Complexities regionsB = new Complexities(chunkB, "foo", 20, 40, 3, 15, template(30), true, null);
    regionsB.globalIntegrity();
    assertNotNull(regionsA.danglingEnd());
    assertEquals(18, regionsA.danglingEnd().getStart());
    assertEquals(19, regionsA.danglingEnd().getEnd());
    assertEquals(ComplexRegion.RegionType.INTERESTING, regionsA.danglingEnd().type());
    assertNotNull(regionsB.danglingStart());
    assertEquals(20, regionsB.danglingStart().getStart());
    assertEquals(21, regionsB.danglingStart().getEnd());
    assertEquals(ComplexRegion.RegionType.INTERESTING, regionsB.danglingStart().type());
    Complexities.fixDangling(null, regionsA);
    Complexities.fixDangling(regionsA, regionsB);
    Complexities.fixDangling(regionsB, null);
    regionsA.globalIntegrity();
    regionsB.globalIntegrity();
    assertNull(regionsA.danglingStart());
    assertNull(regionsB.danglingStart());
    assertNull(regionsA.danglingEnd());
    assertNull(regionsB.danglingEnd());
    assertEquals(1, regionsA.size());
    assertEquals(1, regionsB.size());
    final Iterator<ComplexRegion> it = regionsA.iterator();
    assertTrue(it.hasNext());
    final ComplexRegion merged = it.next();
    assertEquals(18, merged.getStart());
    assertEquals(21, merged.getEnd());
    assertEquals(ComplexRegion.RegionType.COMPLEX, merged.type());
    assertFalse(it.hasNext());

    final Iterator<ComplexRegion> itB = regionsB.iterator();
    assertTrue(itB.hasNext());
    final ComplexRegion mergedB = itB.next();
    assertEquals(18, mergedB.getStart());
    assertEquals(21, mergedB.getEnd());
    assertEquals(ComplexRegion.RegionType.COMPLEX, mergedB.type());
    assertFalse(itB.hasNext());
  }

  //boundary regions become hyper complex
  public void testFixDangling2() {
    final ArrayList<Variant> chunkA = new ArrayList<>();
    chunkA.add(TestUtils.createVariant(17));
    chunkA.add(TestUtils.createVariant(19));
    final Complexities regionsA = new Complexities(chunkA, "foo", 0, 20, 3, 5, template(30), true, null);
    regionsA.globalIntegrity();
    final ArrayList<Variant> chunkB = new ArrayList<>();
    chunkB.add(TestUtils.createVariant(21, true));
    chunkB.add(TestUtils.createVariant(23, true));
    final Complexities regionsB = new Complexities(chunkB, "foo", 20, 40, 3, 5, template(30), true, null);
    regionsB.globalIntegrity();
    assertNotNull(regionsA.danglingEnd());
    assertEquals(17, regionsA.danglingEnd().getStart());
    assertEquals(20, regionsA.danglingEnd().getEnd());
    assertEquals(ComplexRegion.RegionType.COMPLEX, regionsA.danglingEnd().type());
    assertNotNull(regionsB.danglingStart());
    assertEquals(21, regionsB.danglingStart().getStart());
    assertEquals(24, regionsB.danglingStart().getEnd());
    assertEquals(ComplexRegion.RegionType.COMPLEX, regionsB.danglingStart().type());
    Complexities.fixDangling(null, regionsA);
    Complexities.fixDangling(regionsA, regionsB);
    Complexities.fixDangling(regionsB, null);
    regionsA.globalIntegrity();
    regionsB.globalIntegrity();
    assertNull(regionsA.danglingStart());
    assertNotNull(regionsB.danglingStart());
    assertNotNull(regionsA.danglingEnd());
    assertNull(regionsB.danglingEnd());
    assertEquals(1, regionsA.size());
    assertEquals(1, regionsB.size());
    final Iterator<ComplexRegion> itA = regionsA.iterator();
    assertTrue(itA.hasNext());
    final ComplexRegion mergedA = itA.next();
    assertTrue(mergedA == regionsA.danglingEnd());
    assertEquals(17, mergedA.getStart());
    assertEquals(20, mergedA.getEnd());
    assertEquals(ComplexRegion.RegionType.HYPER, mergedA.type());
    assertFalse(itA.hasNext());

    final Iterator<ComplexRegion> itB = regionsB.iterator();
    assertTrue(itB.hasNext());
    final ComplexRegion mergedB = itB.next();
    assertTrue(mergedB == regionsB.danglingStart());
    assertEquals(20, mergedB.getStart());
    assertEquals(24, mergedB.getEnd());
    assertEquals(ComplexRegion.RegionType.HYPER, mergedB.type());
    assertFalse(itB.hasNext());
  }

  //boundary regions extend hyper complex from left
  public void testFixDangling3() {
    final ArrayList<Variant> chunkA = new ArrayList<>();
    chunkA.add(TestUtils.createVariant(13));
    chunkA.add(TestUtils.createVariant(15));
    chunkA.add(TestUtils.createVariant(17));
    chunkA.add(TestUtils.createVariant(19));
    final Complexities regionsA = new Complexities(chunkA, "foo", 0, 20, 3, 5, template(30), true, null);
    regionsA.globalIntegrity();
    final ArrayList<Variant> chunkB = new ArrayList<>();
    chunkB.add(TestUtils.createVariant(21, true));
    chunkB.add(TestUtils.createVariant(23, true));
    final Complexities regionsB = new Complexities(chunkB, "foo", 20, 40, 3, 5, template(30), true, null);
    regionsB.globalIntegrity();
    assertNotNull(regionsA.danglingEnd());
    assertEquals(13, regionsA.danglingEnd().getStart());
    assertEquals(20, regionsA.danglingEnd().getEnd());
    assertEquals(ComplexRegion.RegionType.HYPER, regionsA.danglingEnd().type());
    assertNotNull(regionsB.danglingStart());
    assertEquals(21, regionsB.danglingStart().getStart());
    assertEquals(24, regionsB.danglingStart().getEnd());
    assertEquals(ComplexRegion.RegionType.COMPLEX, regionsB.danglingStart().type());
    Complexities.fixDangling(null, regionsA);
    Complexities.fixDangling(regionsA, regionsB);
    Complexities.fixDangling(regionsB, null);
    regionsA.globalIntegrity();
    regionsB.globalIntegrity();
    assertNull(regionsA.danglingStart());
    assertNotNull(regionsB.danglingStart());
    assertNotNull(regionsA.danglingEnd());
    assertNull(regionsB.danglingEnd());
    assertEquals(1, regionsA.size());
    assertEquals(1, regionsB.size());
    final Iterator<ComplexRegion> itA = regionsA.iterator();
    assertTrue(itA.hasNext());
    final ComplexRegion mergedA = itA.next();
    assertTrue(mergedA == regionsA.danglingEnd());
    assertEquals(13, mergedA.getStart());
    assertEquals(20, mergedA.getEnd());
    assertEquals(ComplexRegion.RegionType.HYPER, mergedA.type());
    assertFalse(itA.hasNext());

    final Iterator<ComplexRegion> itB = regionsB.iterator();
    assertTrue(itB.hasNext());
    final ComplexRegion mergedB = itB.next();
    assertTrue(mergedB == regionsB.danglingStart());
    assertEquals(20, mergedB.getStart());
    assertEquals(24, mergedB.getEnd());
    assertEquals(ComplexRegion.RegionType.HYPER, mergedB.type());
    assertFalse(itB.hasNext());
  }


  //boundary regions extend hyper complex from right
  public void testFixDangling4() {
    final ArrayList<Variant> chunkA = new ArrayList<>();
    chunkA.add(TestUtils.createVariant(17));
    chunkA.add(TestUtils.createVariant(19));
    final Complexities regionsA = new Complexities(chunkA, "foo", 0, 20, 3, 5, template(30), true, null);
    regionsA.globalIntegrity();
    final ArrayList<Variant> chunkB = new ArrayList<>();
    chunkB.add(TestUtils.createVariant(21, true));
    chunkB.add(TestUtils.createVariant(23, true));
    chunkB.add(TestUtils.createVariant(25, true));
    chunkB.add(TestUtils.createVariant(27, true));
    final Complexities regionsB = new Complexities(chunkB, "foo", 20, 40, 3, 5, template(30), true, null);
    regionsB.globalIntegrity();
    assertNotNull(regionsA.danglingEnd());
    assertEquals(17, regionsA.danglingEnd().getStart());
    assertEquals(20, regionsA.danglingEnd().getEnd());
    assertEquals(ComplexRegion.RegionType.COMPLEX, regionsA.danglingEnd().type());
    assertNotNull(regionsB.danglingStart());
    assertEquals(21, regionsB.danglingStart().getStart());
    assertEquals(28, regionsB.danglingStart().getEnd());
    assertEquals(ComplexRegion.RegionType.HYPER, regionsB.danglingStart().type());
    Complexities.fixDangling(null, regionsA);
    Complexities.fixDangling(regionsA, regionsB);
    Complexities.fixDangling(regionsB, null);
    regionsA.globalIntegrity();
    regionsB.globalIntegrity();
    assertNull(regionsA.danglingStart());
    assertNotNull(regionsB.danglingStart());
    assertNotNull(regionsA.danglingEnd());
    assertNull(regionsB.danglingEnd());
    assertEquals(1, regionsA.size());
    assertEquals(1, regionsB.size());
    final Iterator<ComplexRegion> itA = regionsA.iterator();
    assertTrue(itA.hasNext());
    final ComplexRegion mergedA = itA.next();
    assertTrue(mergedA == regionsA.danglingEnd());
    assertEquals(17, mergedA.getStart());
    assertEquals(20, mergedA.getEnd());
    assertEquals(ComplexRegion.RegionType.HYPER, mergedA.type());
    assertFalse(itA.hasNext());

    final Iterator<ComplexRegion> itB = regionsB.iterator();
    assertTrue(itB.hasNext());
    final ComplexRegion mergedB = itB.next();
    assertTrue(mergedB == regionsB.danglingStart());
    assertEquals(20, mergedB.getStart());
    assertEquals(28, mergedB.getEnd());
    assertEquals(ComplexRegion.RegionType.HYPER, mergedB.type());
    assertFalse(itB.hasNext());
  }

  //boundary regions extend hyper complex region across multiple merges
  public void testFixDanglingVeryLongHyperComplex() {
    final ArrayList<Variant> chunkA = new ArrayList<>();
    chunkA.add(TestUtils.createVariant(13));
    chunkA.add(TestUtils.createVariant(15));
    chunkA.add(TestUtils.createVariant(17));
    chunkA.add(TestUtils.createVariant(19));
    final Complexities regionsA = new Complexities(chunkA, "foo", 0, 20, 3, 5, template(30), true, null);
    regionsA.globalIntegrity();
    assertNotNull(regionsA.danglingEnd());
    assertEquals(13, regionsA.danglingEnd().getStart());
    assertEquals(20, regionsA.danglingEnd().getEnd());
    assertEquals(ComplexRegion.RegionType.HYPER, regionsA.danglingEnd().type());

    final ArrayList<Variant> chunkB = new ArrayList<>();
    for (int i = 21; i < 40; i += 2) {
      chunkB.add(TestUtils.createVariant(i, true));
    }
    final Complexities regionsB = new Complexities(chunkB, "foo", 20, 40, 3, 5, template(30), true, null);
    regionsB.globalIntegrity();
    assertNotNull(regionsB.danglingStart());
    assertEquals(21, regionsB.danglingStart().getStart());
    assertEquals(40, regionsB.danglingStart().getEnd());
    assertEquals(ComplexRegion.RegionType.HYPER, regionsB.danglingStart().type());

    final ArrayList<Variant> chunkC = new ArrayList<>();
    chunkC.add(TestUtils.createVariant(41));
    chunkC.add(TestUtils.createVariant(43));
    chunkC.add(TestUtils.createVariant(45));
    chunkC.add(TestUtils.createVariant(47));
    final Complexities regionsC = new Complexities(chunkC, "foo", 40, 60, 3, 5, template(30), true, null);
    regionsC.globalIntegrity();
    assertNotNull(regionsC.danglingStart());
    assertEquals(41, regionsC.danglingStart().getStart());
    assertEquals(48, regionsC.danglingStart().getEnd());
    assertEquals(ComplexRegion.RegionType.HYPER, regionsC.danglingStart().type());

    //merge a and b
    Complexities.fixDangling(null, regionsA);
    Complexities.fixDangling(regionsA, regionsB);
    regionsA.globalIntegrity();
    regionsB.globalIntegrity();
    assertNull(regionsA.danglingStart());
    assertNotNull(regionsB.danglingStart());
    assertNotNull(regionsA.danglingEnd());
    assertNotNull(regionsB.danglingEnd());
    assertEquals(1, regionsA.size());
    assertEquals(1, regionsB.size());
    final Iterator<ComplexRegion> itA = regionsA.iterator();
    assertTrue(itA.hasNext());
    final ComplexRegion mergedA = itA.next();
    assertTrue(mergedA == regionsA.danglingEnd());
    assertEquals(13, mergedA.getStart());
    assertEquals(20, mergedA.getEnd());
    assertEquals(ComplexRegion.RegionType.HYPER, mergedA.type());
    assertFalse(itA.hasNext());

    final Iterator<ComplexRegion> itB = regionsB.iteratorInternal(); //looking at itermediate state
    assertTrue(itB.hasNext());
    final ComplexRegion mergedB = itB.next();
    assertTrue(mergedB == regionsB.danglingStart());
    assertTrue(mergedB == regionsB.danglingEnd());
    assertEquals(20, mergedB.getStart());
    assertEquals(40, mergedB.getEnd());
    assertEquals(ComplexRegion.RegionType.HYPER, mergedB.type());
    assertFalse(itB.hasNext());

    //merge b and c
    Complexities.fixDangling(regionsB, regionsC);
    Complexities.fixDangling(regionsC, null);
    regionsB.globalIntegrity();
    regionsC.globalIntegrity();
    assertNotNull(regionsB.danglingStart());
    assertNotNull(regionsC.danglingStart());
    assertNotNull(regionsB.danglingEnd());
    assertNull(regionsC.danglingEnd());
    assertEquals(1, regionsB.size());
    assertEquals(1, regionsC.size());
    final Iterator<ComplexRegion> itB2 = regionsB.iterator();
    assertTrue(itB2.hasNext());
    final ComplexRegion mergedB2 = itB2.next();
    assertTrue(mergedB2 == regionsB.danglingEnd());
    assertTrue(mergedB2 == regionsB.danglingStart()); //<-- its all for this line really
    assertEquals(20, mergedB2.getStart());
    assertEquals(40, mergedB2.getEnd());
    assertEquals(ComplexRegion.RegionType.HYPER, mergedB2.type());
    assertFalse(itB2.hasNext());

    final Iterator<ComplexRegion> itC = regionsC.iterator();
    assertTrue(itC.hasNext());
    final ComplexRegion mergedC = itC.next();
    assertTrue(mergedC == regionsC.danglingStart());
    assertEquals(40, mergedC.getStart());
    assertEquals(48, mergedC.getEnd());
    assertEquals(ComplexRegion.RegionType.HYPER, mergedC.type());
    assertFalse(itC.hasNext());
  }


  //simple complex calls on left boundary but none on right
  public void testFixDanglingLeftOnly() {
    final ArrayList<Variant> chunkA = new ArrayList<>();
    chunkA.add(TestUtils.createVariant(19, true));
    final Complexities regionsA = new Complexities(chunkA, "foo", 0, 20, 3, 15, template(30), true, null);
    regionsA.globalIntegrity();
    final ArrayList<Variant> chunkB = new ArrayList<>();
    final Complexities regionsB = new Complexities(chunkB, "foo", 20, 40, 3, 15, template(30), true, null);
    regionsB.globalIntegrity();
    assertNotNull(regionsA.danglingEnd());
    assertEquals(19, regionsA.danglingEnd().getStart());
    assertEquals(20, regionsA.danglingEnd().getEnd());
    assertEquals(ComplexRegion.RegionType.COMPLEX, regionsA.danglingEnd().type());
    assertNull(regionsB.danglingStart());
    Complexities.fixDangling(null, regionsA);
    Complexities.fixDangling(regionsA, regionsB);
    Complexities.fixDangling(regionsB, null);
    regionsA.globalIntegrity();
    regionsB.globalIntegrity();
    assertNull(regionsA.danglingStart());
    assertNull(regionsB.danglingStart());
    assertNull(regionsA.danglingEnd());
    assertNull(regionsB.danglingEnd());
    assertEquals(1, regionsA.size());
    assertEquals(0, regionsB.size());
    final Iterator<ComplexRegion> it = regionsA.iterator();
    assertTrue(it.hasNext());
    final ComplexRegion merged = it.next();
    assertEquals(19, merged.getStart());
    assertEquals(20, merged.getEnd());
    assertEquals(ComplexRegion.RegionType.COMPLEX, merged.type());
    assertFalse(it.hasNext());
  }

  //simple complex calls on left but no right region
  public void testFixDanglingLeftRegionOnly() {
    final ArrayList<Variant> chunkA = new ArrayList<>();
    chunkA.add(TestUtils.createVariant(19, true));
    final Complexities regionsA = new Complexities(chunkA, "foo", 0, 20, 3, 15, template(30), true, null);
    regionsA.globalIntegrity();
    assertNotNull(regionsA.danglingEnd());
    assertEquals(19, regionsA.danglingEnd().getStart());
    assertEquals(20, regionsA.danglingEnd().getEnd());
    assertEquals(ComplexRegion.RegionType.COMPLEX, regionsA.danglingEnd().type());
    Complexities.fixDangling(null, regionsA);
    Complexities.fixDangling(regionsA, null);
    regionsA.globalIntegrity();
    assertNull(regionsA.danglingStart());
    assertNull(regionsA.danglingEnd());
    assertEquals(1, regionsA.size());
    final Iterator<ComplexRegion> it = regionsA.iterator();
    assertTrue(it.hasNext());
    final ComplexRegion merged = it.next();
    assertEquals(19, merged.getStart());
    assertEquals(20, merged.getEnd());
    assertEquals(ComplexRegion.RegionType.COMPLEX, merged.type());
    assertFalse(it.hasNext());
  }

  //simple complex calls on right boundary but none on left
  public void testFixDanglingRightOnly() {
    final ArrayList<Variant> chunkA = new ArrayList<>();
    final Complexities regionsA = new Complexities(chunkA, "foo", 0, 20, 3, 15, template(30), true, null);
    regionsA.globalIntegrity();
    final ArrayList<Variant> chunkB = new ArrayList<>();
    chunkB.add(TestUtils.createVariant(21, true));
    final Complexities regionsB = new Complexities(chunkB, "foo", 20, 40, 3, 15, template(30), true, null);
    regionsB.globalIntegrity();
    assertNull(regionsA.danglingEnd());
    assertNotNull(regionsB.danglingStart());
    assertEquals(21, regionsB.danglingStart().getStart());
    assertEquals(22, regionsB.danglingStart().getEnd());
    assertEquals(ComplexRegion.RegionType.COMPLEX, regionsB.danglingStart().type());
    Complexities.fixDangling(null, regionsA);
    Complexities.fixDangling(regionsA, regionsB);
    Complexities.fixDangling(regionsB, null);
    regionsA.globalIntegrity();
    regionsB.globalIntegrity();
    assertNull(regionsA.danglingStart());
    assertNull(regionsB.danglingStart());
    assertNull(regionsA.danglingEnd());
    assertNull(regionsB.danglingEnd());
    assertEquals(0, regionsA.size());
    assertEquals(1, regionsB.size());
    final Iterator<ComplexRegion> it = regionsB.iterator();
    assertTrue(it.hasNext());
    final ComplexRegion merged = it.next();
    assertEquals(21, merged.getStart());
    assertEquals(22, merged.getEnd());
    assertEquals(ComplexRegion.RegionType.COMPLEX, merged.type());
    assertFalse(it.hasNext());
  }

  //simple complex calls on right but no region on left
  public void testFixDanglingRightRegionOnly() {
    final ArrayList<Variant> chunkB = new ArrayList<>();
    chunkB.add(TestUtils.createVariant(21, true));
    final Complexities regionsB = new Complexities(chunkB, "foo", 20, 40, 3, 15, template(30), true, null);
    regionsB.globalIntegrity();
    assertNotNull(regionsB.danglingStart());
    assertEquals(21, regionsB.danglingStart().getStart());
    assertEquals(22, regionsB.danglingStart().getEnd());
    assertEquals(ComplexRegion.RegionType.COMPLEX, regionsB.danglingStart().type());
    Complexities.fixDangling(null, regionsB);
    Complexities.fixDangling(regionsB, null);
    regionsB.globalIntegrity();
    assertNull(regionsB.danglingStart());
    assertNull(regionsB.danglingEnd());
    assertEquals(1, regionsB.size());
    final Iterator<ComplexRegion> it = regionsB.iterator();
    assertTrue(it.hasNext());
    final ComplexRegion merged = it.next();
    assertEquals(21, merged.getStart());
    assertEquals(22, merged.getEnd());
    assertEquals(ComplexRegion.RegionType.COMPLEX, merged.type());
    assertFalse(it.hasNext());
  }

  public void testHyperComplexBoundaryCase() {
    final Complexities regions = new Complexities(Collections.<Variant>emptyList(), "foo", 20, 40, 3, 15, template(30), true, null);
    final LinkedList<ComplexRegion> crs = new LinkedList<>();
    regions.addRegion(crs, 0, 15, true);
    assertEquals(1, crs.size());
    assertEquals(ComplexRegion.RegionType.COMPLEX, crs.peekFirst().type());
  }

  public void testHyperComplexBoundaryCase2() {
    final Complexities regions = new Complexities(Collections.<Variant>emptyList(), "foo", 20, 40, 3, 15, template(30), true, null);
    final LinkedList<ComplexRegion> crs = new LinkedList<>();
    regions.addRegion(crs, 0, 16, true);
    assertEquals(1, crs.size());
    assertEquals(ComplexRegion.RegionType.HYPER, crs.peekFirst().type());
  }

  //simple adjacent complex calls on boundary that aren't close enough
  public void testFixDanglingNoJoin() {
    final ArrayList<Variant> chunkA = new ArrayList<>();
    chunkA.add(TestUtils.createVariant(18, true));
    final Complexities regionsA = new Complexities(chunkA, "foo", 0, 20, 3, 15, template(30), true, null);
    regionsA.globalIntegrity();
    final ArrayList<Variant> chunkB = new ArrayList<>();
    chunkB.add(TestUtils.createVariant(22, true));
    final Complexities regionsB = new Complexities(chunkB, "foo", 20, 40, 3, 15, template(30), true, null);
    regionsB.globalIntegrity();
    assertNotNull(regionsA.danglingEnd());
    assertEquals(18, regionsA.danglingEnd().getStart());
    assertEquals(19, regionsA.danglingEnd().getEnd());
    assertEquals(ComplexRegion.RegionType.COMPLEX, regionsA.danglingEnd().type());
    assertNotNull(regionsB.danglingStart());
    assertEquals(22, regionsB.danglingStart().getStart());
    assertEquals(23, regionsB.danglingStart().getEnd());
    assertEquals(ComplexRegion.RegionType.COMPLEX, regionsB.danglingStart().type());
    Complexities.fixDangling(null, regionsA);
    Complexities.fixDangling(regionsA, regionsB);
    Complexities.fixDangling(regionsB, null);
    regionsA.globalIntegrity();
    regionsB.globalIntegrity();
    assertNull(regionsA.danglingStart());
    assertNull(regionsB.danglingStart());
    assertNull(regionsA.danglingEnd());
    assertNull(regionsB.danglingEnd());
    assertEquals(1, regionsA.size());
    assertEquals(1, regionsB.size());
    final Iterator<ComplexRegion> itA = regionsA.iterator();
    assertTrue(itA.hasNext());
    final ComplexRegion mergedA = itA.next();
    assertEquals(18, mergedA.getStart());
    assertEquals(19, mergedA.getEnd());
    assertEquals(ComplexRegion.RegionType.COMPLEX, mergedA.type());
    assertFalse(itA.hasNext());

    final Iterator<ComplexRegion> itB = regionsB.iterator();
    assertTrue(itB.hasNext());
    final ComplexRegion mergedB = itB.next();
    assertEquals(22, mergedB.getStart());
    assertEquals(23, mergedB.getEnd());
    assertEquals(ComplexRegion.RegionType.COMPLEX, mergedB.type());
    assertFalse(itB.hasNext());

  }

  /*
  This has instead been fixed by making the job scheduler ensure a global ordering of the DANGLING tasks,
  potentially preventing other problems caused by out of order dangling fixing.

  public void testMtIssuesBug1413AndBug1414() throws IOException {
    Diagnostic.setLogStream();
    final SimpleThreadPool stp = new SimpleThreadPool(2, "mytp", true);
    try {
      //change number of iteration > 200,000 for problem to appear faster
      for (int i = 0; i < 20; i++) {
        final ArrayList<Variant> chunkA = new ArrayList<>();
        chunkA.add(TestUtils.createVariant(18));
        final Complexities regionsA = new Complexities(chunkA, "foo", 0, 20, 3, 15, template(30), true, null);
        final ArrayList<Variant> chunkB = new ArrayList<>();
        chunkB.add(TestUtils.createVariant(20));
        chunkB.add(TestUtils.createVariant(22));
        chunkB.add(TestUtils.createVariant(24));
        chunkB.add(TestUtils.createVariant(26));
        chunkB.add(TestUtils.createVariant(28));
        chunkB.add(TestUtils.createVariant(30));
        chunkB.add(TestUtils.createVariant(32));
        chunkB.add(TestUtils.createVariant(34));
        chunkB.add(TestUtils.createVariant(36));
        chunkB.add(TestUtils.createVariant(38));
        final Complexities regionsB = new Complexities(chunkB, "foo", 20, 40, 3, 15, template(30), true, null);
        final ArrayList<Variant> chunkC = new ArrayList<>();
        chunkC.add(TestUtils.createVariant(41));
        final Complexities regionsC = new Complexities(chunkC, "foo", 40, 60, 3, 15, template(30), true, null);

        try {
          stp.execute(new IORunnable() {
            @Override
            public void run() {
              Complexities.fixDangling(regionsA, regionsB);
            }
          });

          stp.execute(new IORunnable() {
            @Override
            public void run() {
              Complexities.fixDangling(regionsB, regionsC);
            }
          });
        } catch (final Exception e) {
          fail(e.getMessage());
        }
        //System.err.println("i= " + i);
      }
    } finally {
      stp.terminate();
    }
  }
*/

  public void testJoin() {
    final ArrayList<Variant> chunkA = new ArrayList<>();
    //                                  0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18
    final byte[] template = {2, 1, 3, 1, 2, 2, 1, 2, 2, 3, 2, 4, 2, 4, 2, 4, 3, 1, 3};
    chunkA.add(TestUtils.createVariant(3));
    chunkA.add(TestUtils.createVariant(8));
    final Complexities regionsA = new Complexities(chunkA, "foo", 0, 20, 4, 15, template, true, null);
    Complexities.fixDangling(null, regionsA);
    Complexities.fixDangling(regionsA, null);
    final Iterator<ComplexRegion> it = regionsA.iterator();
    assertTrue(it.hasNext());
    final ComplexRegion cr = it.next();
    assertEquals(3, cr.getStart());
    assertEquals(9, cr.getEnd());
    assertFalse(it.hasNext());
  }

  public void testMeasureRealWorld() {
    final byte[] template = DNA.stringDNAtoByte("tttccaaaagagtttttcagagaaaaagagatgttgatacatttcaaaaatagga");
    final ArrayList<Variant> chunkA = new ArrayList<>();
    final ArrayList<Variant> chunkB = new ArrayList<>();
    chunkA.add(TestUtils.createVariant(0));
    chunkB.add(TestUtils.createVariant(template.length - 1));
    final Complexities regionsA = new Complexities(chunkA, "foo", 0, 20, 3, 15, template, true, null);
    final Complexities regionsB = new Complexities(chunkB, "foo", 20, 55, 3, 15, template, true, null);
    Complexities.fixDangling(null, regionsA);
    Complexities.fixDangling(regionsA, regionsB);
    Complexities.fixDangling(regionsB, null);
    final Iterator<ComplexRegion> it = regionsA.iterator();
    assertFalse(it.hasNext());
  }
}
