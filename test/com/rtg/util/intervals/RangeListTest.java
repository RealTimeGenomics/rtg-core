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

package com.rtg.util.intervals;

import java.util.ArrayList;
import java.util.List;

import com.rtg.util.intervals.RangeList.RangeData;

import junit.framework.TestCase;

/**
 */
public class RangeListTest extends TestCase {

  public void testRange() {
    RangeList.RangeData<String> r = new RangeData<>(0, 1, "blah");
    assertEquals("blah", r.getMeta().get(0));
    assertTrue(r.isInRange(0));
    assertFalse(r.isInRange(1));
    assertFalse(r.isInRange(-1));
    assertEquals("1-1", r.toString());
    r.addMeta("boo");
    assertEquals(2, r.getMeta().size());
    assertEquals("blah", r.getMeta().get(0));
    assertEquals("boo", r.getMeta().get(1));
    r = new RangeList.RangeData<>(1, 2);
    assertEquals("2-2", r.toString());
    assertEquals(1, r.getStart());
    assertEquals(2, r.getEnd());
  }

  public void testEmptyRange() {
    final RangeList.RangeData<String> r = new RangeList.RangeData<>(0, 0, "blah");
    assertEquals("blah", r.getMeta().get(0));
    assertFalse(r.isInRange(0));
    assertFalse(r.isInRange(1));
    assertFalse(r.isInRange(-1));
    assertEquals("1-0", r.toString());
    r.addMeta("boo");
    assertEquals(2, r.getMeta().size());
    assertEquals("blah", r.getMeta().get(0));
    assertEquals("boo", r.getMeta().get(1));

    final List<RangeList.RangeData<String>> ranges = new ArrayList<>();
    ranges.add(r);
    final RangeList<String> search = new RangeList<>(ranges);

    final List<RangeList.RangeData<String>> newRanges = search.getRangeList();
    assertEquals(0, newRanges.size());

    assertNull(search.find(0));
    assertNull(search.find(1));
    assertNull(search.find(-1));
  }

  public void testRangeSearch() {
    final List<RangeList.RangeData<String>> ranges = new ArrayList<>();
    ranges.add(new RangeList.RangeData<>(1, 4, "1-4"));
    ranges.add(new RangeList.RangeData<>(3, 6, "3-6"));
    ranges.add(new RangeList.RangeData<>(30, 120, "30-120"));
    ranges.add(new RangeData<String>(140, 150));
    final RangeList<String> search = new RangeList<>(ranges);
    List<String> found = search.find(0);
    assertNull(found);
    found = search.find(1);
    assertEquals(1, found.size());
    assertEquals("1-4", found.get(0));
    found = search.find(3);
    assertEquals(2, found.size());
    assertEquals("1-4", found.get(0));
    assertEquals("3-6", found.get(1));
    found = search.find(4);
    assertEquals(1, found.size());
    assertEquals("3-6", found.get(0));
    search.find(6);
    assertNull(null);
    found = search.find(60);
    assertEquals(1, found.size());
    assertEquals("30-120", found.get(0));
    found = search.find(145);
    assertNull(found);
    found = search.find(400);
    assertNull(found);

    found = search.find(Integer.MIN_VALUE);
    assertNull(found);

    found = search.find(Integer.MAX_VALUE - 1);
    assertNull(found);

    found = search.find(Integer.MAX_VALUE);
    assertNull(found);
  }

}


