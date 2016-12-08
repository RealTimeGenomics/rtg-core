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

package com.rtg.variant.realign;

import com.rtg.variant.bayes.complex.ComplexTemplate;

import junit.framework.TestCase;

/**
 */
public class AlignmentEnvironmentGenomeSubstitutionTest extends TestCase {

  public void test() {
    final byte[] template = {1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4};
    final ComplexTemplate cot = new ComplexTemplate(template, "", 11, 13);
    final AlignmentEnvironmentGenomeSubstitution aeg = new AlignmentEnvironmentGenomeSubstitution(3, template.length, cot, new byte[] {1, 3, 2, 1});
    assertEquals(3, aeg.start());
    assertEquals(RealignParamsGenome.SINGLETON.misMatch(), aeg.quality(0));
    assertTrue(!aeg.isInverted());
    assertEquals(template.length + 2 - 3, aeg.subsequenceLength());

    final byte[] exp = {0, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 1, 3, 2, 1, 2, 3, 4, 1, 2, 3, 4, 0};
    for (int i = 0, j = -4; i < exp.length; ++i, ++j) {
      assertEquals(exp[i], aeg.base(j));
    }
  }

  public void testReplaceWithNothing() {
    final byte[] template = {1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4};
    final ComplexTemplate cot = new ComplexTemplate(template, "", 11, 13);
    final AlignmentEnvironmentGenomeSubstitution aeg = new AlignmentEnvironmentGenomeSubstitution(3, template.length, cot, new byte[] {});

    assertEquals(template.length - 2 - 3, aeg.subsequenceLength());

    final byte[] exp = {0, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 2, 3, 4, 1, 2, 3, 4, 0};
    for (int i = 0, j = -4; i < exp.length; ++i, ++j) {
      assertEquals(exp[i], aeg.base(j));
    }
  }

  public void testReplaceWithNothingAtEnd() {
    final byte[] template = {1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4};
    final int l = template.length;
    final ComplexTemplate cot = new ComplexTemplate(template, "", l - 2, l);
    final AlignmentEnvironmentGenomeSubstitution aeg = new AlignmentEnvironmentGenomeSubstitution(3, template.length, cot, new byte[] {});

    assertEquals(template.length - 2 - 3, aeg.subsequenceLength());

    final byte[] exp = {0, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 0};
    for (int i = 0, j = -4; i < exp.length; ++i, ++j) {
      assertEquals(exp[i], aeg.base(j));
    }
  }

  public void testReplaceWithNothingAtStart() {
    final byte[] template = {1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4};
    final ComplexTemplate cot = new ComplexTemplate(template, "", 0, 2);
    final AlignmentEnvironmentGenomeSubstitution aeg = new AlignmentEnvironmentGenomeSubstitution(0, template.length, cot, new byte[] {});

    assertEquals(template.length - 2, aeg.subsequenceLength());

    final byte[] exp = {0, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 0};
    for (int i = 0, j = -1; i < exp.length; ++i, ++j) {
      assertEquals("@" + i, exp[i], aeg.base(j));
    }
  }

  public void testReplaceWithNothingNearStart() {
    final byte[] template = {1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4};
    final ComplexTemplate cot = new ComplexTemplate(template, "", 2, 4);
    final AlignmentEnvironmentGenomeSubstitution aeg = new AlignmentEnvironmentGenomeSubstitution(0, template.length, cot, new byte[] {});

    assertEquals(template.length - 2, aeg.subsequenceLength());

    final byte[] exp = {0, 1, 2, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 0};
    for (int i = 0, j = -1; i < exp.length; ++i, ++j) {
      assertEquals("@" + i, exp[i], aeg.base(j));
    }
  }
  public void testReplaceWithNothingNearStart2() {
    final byte[] template = {1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4};
    final ComplexTemplate cot = new ComplexTemplate(template, "", 2, 4);
    final AlignmentEnvironmentGenomeSubstitution aeg = new AlignmentEnvironmentGenomeSubstitution(3, template.length, cot, new byte[] {});

    assertEquals(template.length - 2 - 3, aeg.subsequenceLength());

    final byte[] exp = {1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 0};
    for (int i = 0, j = -1; i < exp.length; ++i, ++j) {
      assertEquals("@" + i, exp[i], aeg.base(j));
    }
  }
}
