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
