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

package com.rtg.variant.bayes.snp;

import junit.framework.TestCase;

/**
 */
public class DescriptionNoneTest extends TestCase {

  public void test() {
    assertEquals(0, DescriptionNone.SINGLETON.size());
    assertEquals(0, DescriptionNone.SINGLETON.minLength());
    assertEquals(0, DescriptionNone.SINGLETON.maxLength());
    assertFalse(DescriptionNone.SINGLETON.valid(-1));
    assertFalse(DescriptionNone.SINGLETON.valid(0));
  }

}
