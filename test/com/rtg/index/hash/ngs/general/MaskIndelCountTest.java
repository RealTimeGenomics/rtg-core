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
package com.rtg.index.hash.ngs.general;


import java.io.IOException;

import junit.framework.TestCase;

/**
 */
public class MaskIndelCountTest extends TestCase {

  public final void testIndelCount1() throws IOException {
    final Skeleton sk = new Skeleton(5, 3, 2, 1, 1);
    assertEquals(26, MaskIndelCount.indelCount(sk));
  }

  public final void testIndelCount2() throws IOException {
    final Skeleton sk = new Skeleton(5, 3, 2, 0, 1);
    assertEquals(10, MaskIndelCount.indelCount(sk));
  }

  public final void testIndelCount3() throws IOException {
    final Skeleton sk = new Skeleton(1, 1, 0, 0, 1);
    assertEquals(1, MaskIndelCount.indelCount(sk));
  }

  public final void testIndelCount4() throws IOException {
    final Skeleton sk = new Skeleton(6, 3, 3, 0, 1);
    assertEquals(20, MaskIndelCount.indelCount(sk));
  }

  public final void testIndelCount5() throws IOException {
    final Skeleton sk = new Skeleton(6, 3, 3, 1, 1);
    assertEquals(60, MaskIndelCount.indelCount(sk));
  }

  public final void testIndelCount6() throws IOException {
    final Skeleton sk = new Skeleton(6, 3, 3, 1, 2);
    assertEquals(88, MaskIndelCount.indelCount(sk));
  }

}
