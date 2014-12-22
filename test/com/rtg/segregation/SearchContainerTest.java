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

package com.rtg.segregation;

import java.util.HashMap;
import java.util.Map;

import com.rtg.reference.Ploidy;
import com.rtg.reference.ReferenceSequence;
import com.rtg.reference.ReferenceSequenceTest;
import com.rtg.reference.Sex;
import com.rtg.util.Pair;
import com.rtg.util.TestUtils;

import junit.framework.TestCase;

/**
 */
public class SearchContainerTest extends TestCase {


  private static final Map<Pair<Sex, String>, ReferenceSequence> PLOIDYS = new HashMap<>();
  private static final Sex E = Sex.EITHER;

  static {
    PLOIDYS.put(new Pair<>(Sex.EITHER, "foo"), ReferenceSequenceTest.createReferenceSequence(Ploidy.DIPLOID, "1"));
  }

  public void test() throws MismatchingPloidyException {
    final FamilyGt fg0 = FamilyGt.familyPloidy("foo", 10, new String[] {"0/0", "0/1", "0/1", "0/0", "0/0"}, new Sex[] {E, E, E, E, E}, PLOIDYS);
    final SegregationBlock sb0 = new SegregationBlock(fg0);
    final PatternArray pa0 = sb0.patterns();
    final SearchContainer sc0 = new SearchContainer(null, 42.4, null, sb0, pa0, SearchType.New, null, 14);
    sc0.integrity();
    assertEquals(42.4, sc0.score());
    assertTrue(pa0.strictEquals(sc0.pattern()));
    assertTrue(sc0 == sc0.goodContainer());
    assertNull(sc0.previous());
    assertEquals(SearchType.New, sc0.type());
    assertTrue(sb0 == sc0.block());
    assertEquals("14.New<null:42.4 fa: ??? mo: 100", sc0.toString());
    assertEquals(14, sc0.hashCode()); //regression
    assertNull(sc0.xo());

    final FamilyGt fg1 = FamilyGt.familyPloidy("foo", 20, new String[] {"0/1", "0/1", "0/1", "0/0", "1/1"}, new Sex[] {E, E, E, E, E}, PLOIDYS);
    final SegregationBlock sb1 = new SegregationBlock(fg1);
    final PatternArray pa1 = sb1.patterns();
    final SearchContainer sc1 = new SearchContainer(sc0, 42.2, null, sb1, pa1, SearchType.OK, null, 15);
    sc1.integrity();
    assertEquals(42.2, sc1.score());
    assertTrue(pa1.strictEquals(sc1.pattern()));
    assertTrue(sc1 == sc1.goodContainer());
    assertTrue(sc0 == sc1.previous());
    assertEquals(SearchType.OK, sc1.type());
    assertNull(sc1.xo());
    assertTrue(sb1 == sc1.block());
    assertEquals("15.OK<14:42.2 fa: ?01 mo: ?01", sc1.toString());

    final FamilyGt fg2 = FamilyGt.familyPloidy("foo", 20, new String[] {"0/1", "0/1", "0/1", "0/0", "1/1"}, new Sex[] {E, E, E, E, E}, PLOIDYS);
    final SegregationBlock sb2 = new SegregationBlock(fg2);
    final PatternArray pa2 = sb2.patterns();
    final SearchContainer sc2 = new SearchContainer(sc1, 42.0, sc0, sb2, pa2, SearchType.Error, null, 16);
    sc2.integrity();
    assertEquals(42.0, sc2.score());
    assertTrue(pa1.strictEquals(sc2.pattern()));
    assertTrue(sc1 == sc2.previous());
    assertEquals(SearchType.Error, sc2.type());
    assertNull(sc2.xo());
    assertTrue(sb2 == sc2.block());
    assertTrue(sc0 == sc2.goodContainer());
    assertEquals("16.Error<15:42.0 fa: ?01 mo: ?01", sc2.toString());

    assertEquals(0, sc0.compareTo(sc0));
    assertEquals(-1, sc1.compareTo(sc0));
    assertEquals(+1, sc0.compareTo(sc1));

    TestUtils.equalsHashTest(new Object[] {sc0, sc1, sc2});
  }
}
