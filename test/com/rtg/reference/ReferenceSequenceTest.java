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

package com.rtg.reference;

import static com.rtg.util.StringUtils.LS;

import java.util.Iterator;

import com.rtg.util.Pair;
import com.rtg.util.intervals.RegionRestriction;

import junit.framework.TestCase;

/**
 */
public class ReferenceSequenceTest extends TestCase {

  public void test() {
    final ReferenceSequence rs = new ReferenceSequence(true, false, Ploidy.DIPLOID, "xx", null, 0);
    //System.err.println(rs.toString());
    assertFalse(rs.isLinear());
    assertTrue(rs.isSpecified());
    assertEquals(Ploidy.DIPLOID, rs.ploidy());
    assertEquals("xx", rs.name());
    assertFalse(rs.hasDuplicates());
    assertEquals("xx DIPLOID circular 0" + LS, rs.toString());
  }

  public void test1() {
    final ReferenceSequence rs = new ReferenceSequence(false, true, Ploidy.HAPLOID, "yy", "xx", 42);
    final Pair<RegionRestriction, RegionRestriction> d1 = new Pair<>(new RegionRestriction("yy", 0, 15), new RegionRestriction("xx", 7, 99));
    rs.addDuplicate(d1);
    final Pair<RegionRestriction, RegionRestriction> d2 = new Pair<>(new RegionRestriction("xx", 0, 15), new RegionRestriction("yy", 15, 42));
    rs.addDuplicate(d2);
    //System.err.println(rs.toString());
    assertTrue(rs.isLinear());
    assertFalse(rs.isSpecified());
    assertEquals(Ploidy.HAPLOID, rs.ploidy());
    assertEquals(Ploidy.DIPLOID, rs.effectivePloidy(10));
    assertEquals(Ploidy.NONE, rs.effectivePloidy(20));
    assertEquals(Ploidy.HAPLOID, rs.effectivePloidy(43));
    assertEquals("yy", rs.name());
    assertEquals("xx", rs.haploidComplementName());
    assertTrue(rs.hasDuplicates());

    final Iterator<Pair<RegionRestriction, RegionRestriction>> it = rs.duplicates().iterator();
    assertTrue(it.hasNext());
    assertTrue(it.next() == d1);
    assertTrue(it.hasNext());
    assertTrue(it.next() == d2);
    assertFalse(it.hasNext());

    final String exp = ""
      + "yy HAPLOID linear 42 ~=xx" + LS
      + "    yy:1-15  xx:8-99" + LS
      + "    xx:1-15  yy:16-42" + LS
      ;
    assertEquals(exp, rs.toString());
  }

  public void testNonHaploidCheck() {
    final ReferenceSequence rs = new ReferenceSequence(true, true, Ploidy.DIPLOID, "yy", null, 42);
    final Pair<RegionRestriction, RegionRestriction> d1 = new Pair<>(new RegionRestriction("yy", 0, 15), new RegionRestriction("xx", 7, 99));
    try {
      rs.addDuplicate(d1);
      fail();
    } catch (IllegalArgumentException e) {
      assertEquals("Duplicate specified for sequence that isn't haploid.", e.getMessage());
    }
  }

  public static ReferenceSequence createReferenceSequence(Ploidy ploidy, String name) {
    return new ReferenceSequence(true, true, ploidy, name, null, 10000);
  }
}
