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

package com.rtg.variant.bayes.multisample.family;

import com.rtg.variant.bayes.Description;
import com.rtg.variant.bayes.snp.DescriptionCommon;
import com.rtg.variant.bayes.snp.DescriptionSnp;

import junit.framework.TestCase;

/**
 */
public class DescriptionDiseaseTest extends TestCase {

  public void test() {
    final DescriptionCommon desc0 = DescriptionSnp.SINGLETON;
    final Description desc = new DescriptionDisease(desc0);
    assertEquals(5, desc.size());
    assertEquals("NONE", desc.name(0));
    for (int i = 0; i < 4; ++i) {
      assertEquals(desc0.name(i), desc.name(i + 1));
    }
  }

  public void testDifferentLength() {
    final DescriptionCommon desc0 = new DescriptionCommon("X", "YY");
    final Description desc = new DescriptionDisease(desc0);
    assertEquals(3, desc.size());
    assertEquals(1, desc.minLength());
    assertEquals(2, desc.maxLength());
    assertEquals("NONE", desc.name(0));
    for (int i = 0; i < 2; ++i) {
      assertEquals(desc0.name(i), desc.name(i + 1));
    }
  }

}
