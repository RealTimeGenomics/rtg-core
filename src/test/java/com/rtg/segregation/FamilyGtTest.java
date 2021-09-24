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

import junit.framework.TestCase;

/**
 */
public class FamilyGtTest extends TestCase {

  private static final Map<Pair<Sex, String>, ReferenceSequence> PLOIDYS = new HashMap<>();
  private static final Sex E = Sex.EITHER;

  static {
    PLOIDYS.put(new Pair<>(Sex.EITHER, "chr1"), ReferenceSequenceTest.createReferenceSequence(Ploidy.DIPLOID, "1"));
  }

  public void test() throws MismatchingPloidyException {
    final FamilyGt fgt = FamilyGt.familyPloidy("chr1", 42, new String[] {"0/2", "1/2", "0/1", "2/2"}, new Sex[] {E, E, E, E}, PLOIDYS);
    fgt.integrity();
    assertEquals("chr1" + " 42" + " 0_2" + " 1_2" + " 0_1" + " 2_2", fgt.toString());
    assertEquals("chr1", fgt.seq());
    assertEquals(42, fgt.posn());
    assertEquals(2, fgt.length());
    assertEquals("{0/0}", fgt.pattern(0).toString());
    assertEquals("{1/1}", fgt.pattern(1).toString());
    assertTrue(fgt.isMendelian());
    assertFalse(fgt.parentsSingleAllele());
    assertFalse(fgt.isAllHeterozygous());

    assertEquals(0, fgt.father().a());
    assertEquals(2, fgt.father().b());
    assertEquals(1, fgt.mother().a());
    assertEquals(2, fgt.mother().b());
    assertEquals(0, fgt.child(0).a());
    assertEquals(1, fgt.child(0).b());
    assertEquals(2, fgt.child(1).a());
    assertEquals(2, fgt.child(1).b());
  }

  public void testAllIsHeterozygous() throws MismatchingPloidyException {
    checkAllHeterozygous(new String[] {"0/2", "0/2", "0/2", "0/2"}, true);
    checkAllHeterozygous(new String[] {"0/0", "0/0", "0/0", "0/0"}, false);
    checkAllHeterozygous(new String[] {"0/1", "0/2", "0/2", "0/2"}, false);
    checkAllHeterozygous(new String[] {"0/2", "0/1", "0/2", "0/2"}, false);
    checkAllHeterozygous(new String[] {"0/2", "0/2", "0/1", "0/2"}, false);
    checkAllHeterozygous(new String[] {"0/2", "0/2", "0/2", "0/1"}, false);
  }

  public void checkAllHeterozygous(final String[] str, final boolean exp) throws MismatchingPloidyException {
    final FamilyGt fgt = FamilyGt.familyPloidy("chr1", 42, str, new Sex[] {E, E, E, E}, PLOIDYS);
    fgt.integrity();
    assertEquals(exp, fgt.isAllHeterozygous());
  }

  public void testIsHomozygous() throws MismatchingPloidyException {
    final FamilyGt fgt = FamilyGt.familyPloidy("chr1", 42, new String[] {"0/0", "1/1", "0/1", "2/2"}, new Sex[] {E, E, E, E}, PLOIDYS);
    fgt.integrity();
    assertFalse(fgt.isMendelian());
    assertTrue(fgt.parentsSingleAllele());
  }

  public void testNotMendelian() throws MismatchingPloidyException {
    final FamilyGt fgt = FamilyGt.familyPloidy("chr1", 42, new String[] {"0/1", "1/1", "0/1", "2/2"}, new Sex[] {E, E, E, E}, PLOIDYS);
    fgt.integrity();
    assertFalse(fgt.isMendelian());
  }
}
