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
public class LabellingTest extends TestCase {
  private static final Map<Pair<Sex, String>, ReferenceSequence> PLOIDYS = new HashMap<>();
  private static final Sex E = Sex.EITHER;

  static {
    PLOIDYS.put(new Pair<>(Sex.EITHER, "chr1"), ReferenceSequenceTest.createReferenceSequence(Ploidy.DIPLOID, "1"));
  }


  public void testEquivalent() {
    checkEquivalent(0, 0, 0, 0, true);
    checkEquivalent(1, 0, 1, 0, true);
    checkEquivalent(1, 0, 0, 1, true);
    checkEquivalent(0, 1, 1, 0, true);
    checkEquivalent(0, 1, 0, 1, true);

    checkEquivalent(0, 0, 0, 1, false);
    checkEquivalent(0, 0, 1, 0, false);
    checkEquivalent(0, 1, 0, 0, false);
    checkEquivalent(1, 0, 0, 0, false);
  }

  private void checkEquivalent(final int a, final int b, final int c, final int d, final boolean exp) {
    assertEquals(exp, Labelling.equivalent(new GType(a, b, null), new GType(c, d, null)));
  }

  public void testASame() {
    checkASame(0, 0, 0, 0, true);
    checkASame(1, 0, 1, 0, true);
    checkASame(1, 0, 0, 1, true);
    checkASame(0, 1, 1, 0, true);
    checkASame(0, 1, 0, 1, true);

    checkASame(0, 0, 0, 1, true);
    checkASame(0, 0, 1, 0, true);
    checkASame(0, 1, 0, 0, true);
    checkASame(1, 0, 0, 0, false);

    checkASame(0, 0, 1, 2, false);

  }

  private void checkASame(final int a, final int b, final int c, final int d, final boolean exp) {
    assertEquals(exp, Labelling.aSame(new GType(a, b, null), new GType(c, d, null)));
  }

  public void testLabelPhasing() throws MismatchingPloidyException {
    final FamilyGt fgt = FamilyGt.familyPloidy("chr1", 42, new String[] {"0/1", "0/1", "0/0", "0/1", "0/1", "1/1"}, new Sex[] {E, E, E, E, E, E}, PLOIDYS);
    //System.err.println(fgt.pattern().toString());
    final PatternArray pattern = new PatternArray("0101", "0011");
    final Labelling lab = new Labelling(fgt, pattern);
    lab.globalIntegrity();

    check(lab.labelPhasing(1), "1_0");
    check(lab.labelPhasing(2), "0_1");

    check(lab.phaseChild(1), "1_0");
    check(lab.phaseChild(2), "0_1");

    check(lab.phaseFather(), "0_1");
    check(lab.phaseMother(), "0_1");
  }

  public void testLabelPhasingFaInvalid() throws MismatchingPloidyException {
    final FamilyGt fgt = FamilyGt.familyPloidy("chr1", 42, new String[] {"0/1", "0/1", "0/0", "0/1", "0/1", "1/1"}, new Sex[] {E, E, E, E, E, E}, PLOIDYS);
    //System.err.println(fgt.pattern().toString());
    final PatternArray pattern = new PatternArray("?101", "0011");
    final Labelling lab = new Labelling(fgt, pattern);
    lab.globalIntegrity();

    check(lab.labelPhasing(1), "1_0");
    check(lab.labelPhasing(2), "0_1");

    check(lab.phaseChild(1), "1_0");
    check(lab.phaseChild(2), "0_1");

    check(lab.phaseFather(), "0_1");
    check(lab.phaseMother(), "0_1");
  }

  public void testLabelPhasingMoInvalidX() throws MismatchingPloidyException {
    final FamilyGt fgt = FamilyGt.familyPloidy("chr1", 42, new String[] {"0/1", "0/1", "0/0", "0/1", "0/1", "1/1"}, new Sex[] {E, E, E, E, E, E}, PLOIDYS);
    //System.err.println(fgt.pattern().toString());
    final PatternArray pattern = new PatternArray("0101", "?100");
    final Labelling lab = new Labelling(fgt, pattern);
    lab.globalIntegrity();

    check(lab.labelPhasing(1), "1_0");
    check(lab.labelPhasing(2), "0_1");

    check(lab.phaseChild(1), "1_0");
    check(lab.phaseChild(2), "0_1");

    check(lab.phaseFather(), "0_1");
    check(lab.phaseMother(), "1_0");
  }

  public void testLabelPhasingTricky() throws MismatchingPloidyException {
    final FamilyGt fgt = FamilyGt.familyPloidy("chr1", 42, new String[] {"0/1", "0/1", "0/0", "0/1"}, new Sex[] {E, E, E, E}, PLOIDYS);
    //System.err.println(fgt.pattern().toString());
    final PatternArray pattern = new PatternArray("01", "00");
    final Labelling lab = new Labelling(fgt, pattern);
    lab.globalIntegrity();

    check(lab.labelPhasing(1), "1_0");

    check(lab.phaseChild(1), "1_0");

    check(lab.phaseFather(), "0_1");
    check(lab.phaseMother(), "0_1");
  }

  public void testLabelPhasingMoInvalid() throws MismatchingPloidyException {
    final FamilyGt fgt = FamilyGt.familyPloidy("chr1", 42, new String[] {"0/1", "0/1", "0/0", "0/1", "0/1", "1/1"}, new Sex[] {E, E, E, E, E, E}, PLOIDYS);
    //System.err.println(fgt.pattern().toString());
    final PatternArray pattern = new PatternArray("0101", "?011");
    final Labelling lab = new Labelling(fgt, pattern);
    lab.globalIntegrity();

    check(lab.labelPhasing(1), "1_0");
    check(lab.labelPhasing(2), "0_1");

    check(lab.phaseChild(1), "1_0");
    check(lab.phaseChild(2), "0_1");

    check(lab.phaseFather(), "0_1");
    check(lab.phaseMother(), "0_1");
  }

  public void testPhaseChild() throws MismatchingPloidyException {
    final FamilyGt fgt = FamilyGt.familyPloidy("chr1", 42, new String[] {"0/2", "1/2", "0/1", "2/2", "1/0", "0/2", "2/0", "1/2", "2/1"}, new Sex[] {E, E, E, E, E, E, E, E, E}, PLOIDYS);
    //System.err.println(fgt.pattern().toString());
    final PatternArray pattern = new PatternArray("0100011", "0101100");
    final Labelling lab = new Labelling(fgt, pattern);
    lab.globalIntegrity();
    check(lab.phaseChild(0), "0_1");
    check(lab.phaseChild(2), "0_1");
    check(lab.phaseChild(3), "0_2");
    check(lab.phaseChild(4), "0_2");
    check(lab.phaseChild(5), "2_1");
    check(lab.phaseChild(6), "2_1");

    check(lab.phaseFather(), "0_2");
    check(lab.phaseMother(), "1_2");
  }

  public void testPhaseChildEachLabelOnce() throws MismatchingPloidyException {
    final FamilyGt fgt = FamilyGt.familyPloidy("chr1", 42, new String[] {"0/2", "1/2", "0/1", "2/2", "1/2"}, new Sex[] {E, E, E, E, E}, PLOIDYS);
    final PatternArray pattern = new PatternArray("011", "010");
    final Labelling lab = new Labelling(fgt, pattern);
    lab.globalIntegrity();
    check(lab.phaseChild(0), "0_1");
    check(lab.phaseChild(2), "2_1");

    check(lab.phaseFather(), "0_2");
    check(lab.phaseMother(), "1_2");
  }

  public void testPhaseChildInvalid() throws MismatchingPloidyException {
    final FamilyGt fgt = FamilyGt.familyPloidy("chr1", 42, new String[] {"0/1", "1/0", "0/0", "0/1"}, new Sex[] {E, E, E, E}, PLOIDYS);
    final PatternArray pattern = new PatternArray("0?", "01");
    final Labelling lab = new Labelling(fgt, pattern);
    lab.globalIntegrity();
    check(lab.phaseChild(1), null);

    check(lab.phaseFather(), "0_1");
    check(lab.phaseMother(), "0_1");
  }

  public void testPhaseChildOneLabelNull() throws MismatchingPloidyException {
    final FamilyGt fgt = FamilyGt.familyPloidy("chr1", 42, new String[] {"0/2", "1/2", "0/1", "2/2", "1/0", "0/2", "2/0"}, new Sex[] {E, E, E, E, E, E, E}, PLOIDYS);
    //System.err.println(fgt.pattern().toString());
    final PatternArray pattern = new PatternArray("01000", "01011");
    final Labelling lab = new Labelling(fgt, pattern);
    lab.globalIntegrity();
    check(lab.phaseChild(0), "0_1");
    check(lab.phaseChild(2), "0_1");
    check(lab.phaseChild(3), "0_2");
    check(lab.phaseChild(4), "0_2");

    check(lab.phaseFather(), "0_2");
    check(lab.phaseMother(), "1_2");
  }

  //Bug seen in real case - should be able to phase and doesnt
  public void testBug() throws MismatchingPloidyException {
    final FamilyGt fgt = FamilyGt.familyPloidy("chr1", 42, new String[] {"0/2", "0/0", "0/0", "0/2", "0/2"}, new Sex[] {E, E, E, E, E}, PLOIDYS);
    //System.err.println(fgt.pattern().toString());
    final PatternArray pattern = new PatternArray("100", "110");
    final Labelling lab = new Labelling(fgt, pattern);
    lab.globalIntegrity();
    check(lab.phaseChild(2), "2_0"); //the buggy case
    check(lab.phaseChild(1), "2_0");

    check(lab.phaseFather(), "2_0");
  }

  private void check(GType p0, String exp) {
    if (exp == null) {
      assertNull(p0);
      return;
    }
    assertEquals(exp, p0.toString());
  }

  public void testIsIn() {
    final GType gt = new GType(1, 2, null);
    assertFalse(Labelling.isIn(0, gt));
    assertTrue(Labelling.isIn(1, gt));
    assertTrue(Labelling.isIn(2, gt));
  }

  public void testSimpleCase() {
    check(1, 0, 0, 0, 1, false);
    check(1, 0, 1, 0, 0, true);
    check(1, 0, 1, 2, 3, true);
    check(0, 0, 1, 0, 0, null);
  }

  private void check(final int x, final int a, final int b, final int c, final int d, final Boolean exp) {
    final GType father = new GType(a, b, null);
    final GType mother = new GType(c, d, null);
    final GType primary = new GType(4, 5, null);
    final GType secondary = new GType(5, 4, null);
    final GType actual = Labelling.simpleCase(x, father, mother, primary, secondary);
    if (exp == null) {
      assertNull(actual);
      return;
    }
    assertEquals(exp.booleanValue(), actual.toString().equals("4_5"));
  }

  public void testSimplePhasing() {
    checkSimple(0, 0, 0, 1, 0, 1, true);
    checkSimple(0, 1, 0, 0, 0, 1, false);
    checkSimple(0, 1, 0, 1, 0, 1, null);
    checkSimple(0, 1, 0, 2, 0, 1, false);
  }

  private void checkSimple(final int a, final int b, final int c, final int d, final int e, final int f, final Boolean exp) {
    final GType father = new GType(a, b, null);
    final GType mother = new GType(c, d, null);
    final GType child = new GType(e, f, null);
    final GType childAlt = new GType(f, e, null);
    final GType actual = Labelling.simplePhasing(father, mother, child);
    if (exp == null) {
      assertNull(actual);
      return;
    }
    assertEquals((exp ? child : childAlt).toString(), actual.toString());
  }


}
