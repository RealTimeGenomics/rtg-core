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

import junit.framework.TestCase;

/**
 */
public class AlignmentEnvironmentGenomeTest extends TestCase {

  public void test() {
    final byte[] template = {1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4};
    final AlignmentEnvironmentGenome aeg = new AlignmentEnvironmentGenome(3, template.length, template);

    assertEquals(template.length - 3, aeg.subsequenceLength());

    final byte[] exp = {0, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 0};
    for (int i = 0, j = -4; i < exp.length; ++i, ++j) {
      assertEquals(exp[i], aeg.base(j));
    }
  }
}
