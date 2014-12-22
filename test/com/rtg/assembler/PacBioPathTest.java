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

package com.rtg.assembler;

import java.util.Arrays;

import junit.framework.TestCase;

/**
 */
public class PacBioPathTest extends TestCase {
  public void testToPath() {
    PacBioPath path = new PacBioPath(null, new PartialAlignment(2, 10, 20, 1, 1, 11));
    assertEquals(Arrays.asList(1L), path.toPath());
    assertEquals(2, path.score());
    assertFalse(path.mIsPrefix);
    path.mIsPrefix = true;
    path = new PacBioPath(path, new PartialAlignment(4, 15, 25, 5, 1, 11));
    assertFalse(path.mIsPrefix);
    assertTrue(path.mPrevious.mIsPrefix);
    assertEquals(Arrays.asList(1L, 5L), path.toPath());
    assertEquals(6, path.score());
    assertEquals(", 1, 5 isDupe=false score=6 readLength=15[10-25]", path.toString());
  }
}
