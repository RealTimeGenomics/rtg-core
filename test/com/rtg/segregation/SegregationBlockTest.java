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
import com.rtg.util.Utils;

import junit.framework.TestCase;

/**
 */
public class SegregationBlockTest extends TestCase {

  private static final Map<Pair<Sex, String>, ReferenceSequence> PLOIDYS = new HashMap<>();
  private static final Sex E = Sex.EITHER;

  static {
    PLOIDYS.put(new Pair<>(Sex.EITHER, "foo"), ReferenceSequenceTest.createReferenceSequence(Ploidy.DIPLOID, "1"));
  }

  public void test() throws MismatchingPloidyException {
    final FamilyGt fgt0 = FamilyGt.familyPloidy("foo", 42, new String[] {"0/0", "0/1", "0/0", "0/1"}, new Sex[] {E, E, E, E}, PLOIDYS);
    fgt0.integrity();
    final SegregationBlock sgb = new SegregationBlock(fgt0);
    sgb.integrity();
    assertEquals(1, sgb.count());
    assertEquals(42, sgb.start());
    assertEquals(42, sgb.end());
    assertEquals("foo 42 - 42 1 fa ?? mo 01", sgb.toString());

    final FamilyGt fgt1 = FamilyGt.familyPloidy("foo", 54, new String[] {"0/1", "0/0", "0/0", "0/1"}, new Sex[] {E, E, E, E}, PLOIDYS);
    assertTrue(sgb.extend(fgt1));
    sgb.integrity();
    assertEquals(2, sgb.count());
    assertEquals(42, sgb.start());
    assertEquals(54, sgb.end());
    assertEquals("foo 42 - 54 2 fa 01 mo 01", sgb.toString());

    final FamilyGt fgt2 = FamilyGt.familyPloidy("foo", 64, new String[] {"0/1", "2/3", "0/2", "0/2"}, new Sex[] {E, E, E, E}, PLOIDYS);
    assertFalse(sgb.extend(fgt2));
    sgb.integrity();
    assertEquals(2, sgb.count());
    assertEquals(42, sgb.start());
    assertEquals(54, sgb.end());
    assertEquals("foo 42 - 54 2 fa 01 mo 01", sgb.toString());
  }

  public void testCrossover() throws MismatchingPloidyException {
    final String[] s0a = {"0/0", "0/1", "0/1", "0/0"};
    final String[] s0b = {"0/1", "0/1", "0/0", "0/0"};
    final String[] s1a = {"0/0", "0/0", "0/1", "0/0"};
    final String[] s1b = {"0/1", "0/1", "0/0", "0/0"};
    final String exp = "1 2 3";
    checkCrossOver(s0a, s0b, s1a, s1b, exp);
  }

  public void testCrossoverFirstChild() throws MismatchingPloidyException {
    final String[] s0a = {"0/0", "0/1", "0/1", "0/0"};
    final String[] s0b = {"0/1", "0/1", "0/0", "0/0"};
    final String[] s1a = {"0/1", "0/1", "0/1", "0/0"};
    final String[] s1b = {"0/1", "0/1", "0/0", "0/0"};
    final String exp = "0 2 1";
    checkCrossOver(s0a, s0b, s1a, s1b, exp);
  }

  public void testCrossoverComplement() throws MismatchingPloidyException {
    final String[] s0a = {"0/0", "0/1", "0/1", "0/0"};
    final String[] s0b = {"0/1", "0/1", "0/0", "0/0"};
    final String[] s1a = {"0/1", "0/0", "0/0", "0/0"};
    final String[] s1b = {"0/1", "0/1", "0/0", "0/0"};
    final String exp = "3 2 1";
    checkCrossOver(s0a, s0b, s1a, s1b, exp);
  }

  public void testCrossoverEqual() throws MismatchingPloidyException {
    final String[] s0a = {"0/0", "0/1", "0/1", "0/0"};
    final String[] s0b = {"0/1", "0/1", "0/0", "0/0"};
    final String[] s1a = {"0/0", "0/1", "0/1", "0/0"};
    final String[] s1b = {"0/1", "0/1", "0/0", "0/0"};
    checkCrossOver(s0a, s0b, s1a, s1b, null);
  }

  public void testCrossoverBothDifferent() throws MismatchingPloidyException {
    final String[] s0a = {"0/0", "0/1", "0/1", "0/0"};
    final String[] s0b = {"0/1", "0/1", "0/0", "0/0"};
    final String[] s1a = {"0/0", "0/0", "0/1", "0/0"};
    final String[] s1b = {"0/1", "0/1", "0/0", "0/1"};
    checkCrossOver(s0a, s0b, s1a, s1b, null);
  }

  private void checkCrossOver(final String[] s0a, final String[] s0b, final String[] s1a, final String[] s1b, final String exp) throws MismatchingPloidyException {
    checkXO(
        Utils.append(new String[] {"0/0", "0/1"}, s0a),
        Utils.append(new String[] {"0/1", "0/0"}, s0b),
        Utils.append(new String[] {"0/0", "0/1"}, s1a),
        Utils.append(new String[] {"0/1", "0/0"}, s1b),
        exp == null ? null : "mo " + exp
        );
    checkXO(
        Utils.append(new String[] {"0/1", "0/0"}, s0a),
        Utils.append(new String[] {"0/0", "0/1"}, s0b),
        Utils.append(new String[] {"0/1", "0/0"}, s1a),
        Utils.append(new String[] {"0/0", "0/1"}, s1b),
        exp == null ? null : "fa " + exp
        );
  }

  private void checkXO(final String[] s0a, final String[] s0b, final String[] s1a, final String[] s1b, final String exp) throws MismatchingPloidyException {
    final String seq = "foo";
    final FamilyGt fgt0a = FamilyGt.familyPloidy(seq, 42, s0a, new Sex[] {E, E, E, E, E, E}, PLOIDYS);
    final FamilyGt fgt0b = FamilyGt.familyPloidy(seq, 54, s0b, new Sex[] {E, E, E, E, E, E}, PLOIDYS);
    final SegregationBlock sgb0 = new SegregationBlock(fgt0a);
    assertTrue(sgb0.extend(fgt0b));
    final FamilyGt fgt1a = FamilyGt.familyPloidy(seq, 64, s1a, new Sex[] {E, E, E, E, E, E}, PLOIDYS);
    final FamilyGt fgt1b = FamilyGt.familyPloidy(seq, 72, s1b, new Sex[] {E, E, E, E, E, E}, PLOIDYS);
    final SegregationBlock sgb1 = new SegregationBlock(fgt1a);
    assertTrue(sgb1.extend(fgt1b));
    //System.err.println(sgb0.toString());
    //System.err.println(sgb1.toString());
    if (exp == null) {
      assertNull(PatternArray.crossover(sgb0.patterns(), sgb1.patterns(), false));
    } else {
      assertEquals(exp, PatternArray.crossover(sgb0.patterns(), sgb1.patterns(), false).toString());
    }
  }
}
