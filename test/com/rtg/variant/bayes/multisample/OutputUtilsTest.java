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
import java.util.Collections;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

import com.rtg.variant.Variant;
import com.rtg.variant.Variant.VariantFilter;

import junit.framework.TestCase;

/**
 */
public class OutputUtilsTest extends TestCase {

  public void testFilterEmpty() {
    final ArrayList<Variant> chunk = new ArrayList<>();
    final Complexities regions = new Complexities(chunk, null, 0, 50, 3, 15, ComplexitiesTest.template(30), true, null);

    Complexities.fixDangling(null, regions);
    Complexities.fixDangling(regions, null);

    final List<Variant> filter = OutputUtils.nonShadowed(chunk, regions);
    assertEquals(0, filter.size());
  }

  public void testFilterSimple() {
    final ArrayList<Variant> chunk = new ArrayList<>();
    chunk.add(TestUtils.createVariant(4));
    chunk.add(TestUtils.createVariant(5));
    chunk.add(TestUtils.createVariant(6));
    chunk.add(TestUtils.createVariant(10));
    chunk.add(TestUtils.createVariant(15));
    final Complexities regions = new Complexities(chunk, "foo", 0, 50, 3, 15, ComplexitiesTest.template(30), false, null);

    Complexities.fixDangling(null, regions);
    Complexities.fixDangling(regions, null);

    final List<Variant> filter = OutputUtils.nonShadowed(chunk, regions);
    assertEquals(2, filter.size());
    final Iterator<Variant> itf = filter.iterator();
    assertTrue(itf.hasNext());
    final Variant call1 = itf.next();
    assertEquals(10, call1.getLocus().getStart());
    assertEquals(11, call1.getLocus().getEnd());
    assertTrue(itf.hasNext());
    final Variant call2 = itf.next();
    assertEquals(15, call2.getLocus().getStart());
    assertEquals(16, call2.getLocus().getEnd());
    assertFalse(itf.hasNext());

  }

  //Two complex regions
  public void testFilterTwoComplex() {
    final ArrayList<Variant> chunk = new ArrayList<>();
    chunk.add(TestUtils.createVariant(4));
    chunk.add(TestUtils.createVariant(5));
    chunk.add(TestUtils.createVariant(6));
    chunk.add(TestUtils.createVariant(10));
    chunk.add(TestUtils.createVariant(15));
    chunk.add(TestUtils.createVariant(17));
    final Complexities regions = new Complexities(chunk, "foo", 0, 50, 3, 15, ComplexitiesTest.template(30), false, null);

    Complexities.fixDangling(null, regions);
    Complexities.fixDangling(regions, null);

    final List<Variant> filter = OutputUtils.nonShadowed(chunk, regions);
    assertEquals(1, filter.size());
    final Iterator<Variant> itf = filter.iterator();
    assertTrue(itf.hasNext());
    final Variant call1 = itf.next();
    assertEquals(10, call1.getLocus().getStart());
    assertEquals(11, call1.getLocus().getEnd());
    assertFalse(itf.hasNext());

  }

  //Two hyper-complex regions
  public void testFilterTwoHyper() {
    final ArrayList<Variant> chunk = new ArrayList<>();
    chunk.add(TestUtils.createVariant(4));
    chunk.add(TestUtils.createVariant(5));
    chunk.add(TestUtils.createVariant(7));
    chunk.add(TestUtils.createVariant(11));
    chunk.add(TestUtils.createVariant(15));
    chunk.add(TestUtils.createVariant(17));
    chunk.add(TestUtils.createVariant(19));
    final Complexities regions = new Complexities(chunk, "foo", 0, 50, 3, 3, ComplexitiesTest.template(30), true, null);

    Complexities.fixDangling(null, regions);
    Complexities.fixDangling(regions, null);

    final List<Variant> filter = OutputUtils.nonShadowed(chunk, regions);
    assertEquals(7, filter.size());
    final Iterator<Variant> itf = filter.iterator();
    assertTrue(itf.hasNext());
    check(itf, 4, 5, true);
    check(itf, 5, 6, true);
    check(itf, 7, 8, true);
    check(itf, 11, 12, false);
    check(itf, 15, 16, true);
    check(itf, 17, 18, true);
    check(itf, 19, 20, true);

    assertFalse(itf.hasNext());

  }

  //boundary regions extend hyper complex region across multiple merges
  public void testFixDanglingVeryLongHyperComplex() {
    final ArrayList<Variant> chunkA = new ArrayList<>();
    chunkA.add(TestUtils.createVariant(13));
    chunkA.add(TestUtils.createVariant(15));
    chunkA.add(TestUtils.createVariant(17));
    chunkA.add(TestUtils.createVariant(19));
    final Complexities regionsA = new Complexities(chunkA, "foo", 0, 20, 3, 5, ComplexitiesTest.template(30), true, null);

    final ArrayList<Variant> chunkB = new ArrayList<>();
    for (int i = 21; i < 40; i += 2) {
      chunkB.add(TestUtils.createVariant(i, true));
    }
    final Complexities regionsB = new Complexities(chunkB, "foo", 20, 40, 3, 5, ComplexitiesTest.template(30), true, null);

    final ArrayList<Variant> chunkC = new ArrayList<>();
    chunkC.add(TestUtils.createVariant(41));
    chunkC.add(TestUtils.createVariant(43));
    chunkC.add(TestUtils.createVariant(45));
    chunkC.add(TestUtils.createVariant(47));
    final Complexities regionsC = new Complexities(chunkC, "foo", 40, 60, 3, 5, ComplexitiesTest.template(30), true, null);

    //merge a and b and c
    Complexities.fixDangling(null, regionsA);
    Complexities.fixDangling(regionsA, regionsB);
    Complexities.fixDangling(regionsB, regionsC);
    Complexities.fixDangling(regionsC, null);

    final List<Variant> filterA = OutputUtils.nonShadowed(chunkA, regionsA);
    assertEquals(4, filterA.size());
    final Iterator<Variant> itA = filterA.iterator();
    check(itA, 13, 14, true);
    check(itA, 15, 16, true);
    check(itA, 17, 18, true);
    check(itA, 19, 20, true);
    assertFalse(itA.hasNext());

    final List<Variant> filterB = OutputUtils.nonShadowed(chunkB, regionsB);
    assertEquals(0, filterB.size());

    final List<Variant> filterC = OutputUtils.nonShadowed(chunkC, regionsC);
    assertEquals(4, filterC.size());
    final Iterator<Variant> itC = filterC.iterator();
    check(itC, 41, 42, true);
    check(itC, 43, 44, true);
    check(itC, 45, 46, true);
    check(itC, 47, 48, true);
    assertFalse(itC.hasNext());


  }

  void check(Iterator<Variant> it, int start, int end, boolean hyper) {
    final Variant call = it.next();
    assertEquals(start, call.getLocus().getStart());
    assertEquals(end, call.getLocus().getEnd());
    assertEquals(hyper, call.isFiltered(VariantFilter.HYPER_COMPLEX));
  }

  public void testMergeEmpty() {
    final List<Variant> merge = OutputUtils.merge(Collections.<Variant>emptyList(), Collections.<Variant>emptyList());
    assertEquals(0, merge.size());
  }

  public void testMergeEmptyLeft() {
    final List<Variant> right = new LinkedList<>();
    right.add(TestUtils.createVariant(0));
    right.add(TestUtils.createVariant(3));

    final List<Variant> merge = OutputUtils.merge(Collections.<Variant>emptyList(), right);
    assertEquals(2, merge.size());
    final Iterator<Variant> it = merge.iterator();
    assertEquals(0, it.next().getLocus().getStart());
    assertEquals(3, it.next().getLocus().getStart());
    assertFalse(it.hasNext());
  }

  public void testMergeEmptyRight() {
    final List<Variant> left = new LinkedList<>();
    left.add(TestUtils.createVariant(0));
    left.add(TestUtils.createVariant(3));

    final List<Variant> merge = OutputUtils.merge(left, Collections.<Variant>emptyList());
    assertEquals(2, merge.size());
    final Iterator<Variant> it = merge.iterator();
    assertEquals(0, it.next().getLocus().getStart());
    assertEquals(3, it.next().getLocus().getStart());
    assertFalse(it.hasNext());
  }

  public void testMerge() {
    final List<Variant> left = new LinkedList<>();
    left.add(TestUtils.createVariant(0));
    left.add(TestUtils.createVariant(3));
    final List<Variant> right = new LinkedList<>();
    right.add(TestUtils.createVariant(2));
    right.add(TestUtils.createVariant(5));

    final List<Variant> merge = OutputUtils.merge(left, right);
    assertEquals(4, merge.size());
    final Iterator<Variant> it = merge.iterator();
    assertEquals(0, it.next().getLocus().getStart());
    assertEquals(2, it.next().getLocus().getStart());
    assertEquals(3, it.next().getLocus().getStart());
    assertEquals(5, it.next().getLocus().getStart());
    assertFalse(it.hasNext());
  }

  public void testMergeSafeOverlap() {
    final List<Variant> left = new LinkedList<>();
    left.add(TestUtils.createVariant(3));
    final List<Variant> right = new LinkedList<>();
    right.add(TestUtils.createVariant(3, 3, false, 0));
    final List<Variant> merge = OutputUtils.merge(left, right);
    assertEquals(2, merge.size());
    final Iterator<Variant> it = merge.iterator();
    assertEquals(3, it.next().getLocus().getStart());
    assertEquals(3, it.next().getLocus().getStart());
    assertFalse(it.hasNext());
  }

  public void testMergeBadOverlap() {
    final List<Variant> left = new LinkedList<>();
    left.add(TestUtils.createVariant(3, 4, false, 0));
    final List<Variant> right = new LinkedList<>();
    right.add(TestUtils.createVariant(3, 5, false, 0));
    final List<Variant> merge = OutputUtils.merge(left, right);
    assertEquals(2, merge.size()); // We now support merging multiple variants with the same start position
    final Iterator<Variant> it = merge.iterator();
    assertEquals(4, it.next().getLocus().getEnd());
    assertEquals(5, it.next().getLocus().getEnd());
    assertFalse(it.hasNext());
  }

}
