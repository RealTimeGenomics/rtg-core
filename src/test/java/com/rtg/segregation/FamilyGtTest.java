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
