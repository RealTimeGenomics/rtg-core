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

package com.rtg.reader;

import junit.framework.TestCase;

/**
 */
public class BestSumReadTrimmerTest extends TestCase {

  public void testVeryGoodQual() {
    final byte[] quals = {73, 73, 73, 73, 73, 73, 73, 73, 73};
    final BestSumReadTrimmer bsrt = new BestSumReadTrimmer(15);

    assertEquals(quals.length, bsrt.getTrimPosition(quals, quals.length));
    assertEquals(8, bsrt.getTrimPosition(quals, 8));
  }

  public void testHalfGoodQual() {
    final byte[] quals = {73, 73, 73, 73, 13, 13, 13, 13};
    BestSumReadTrimmer bsrt = new BestSumReadTrimmer(15);

    assertEquals(4, bsrt.getTrimPosition(quals, quals.length));
    assertEquals(4, bsrt.getTrimPosition(quals, 6));
    bsrt = new BestSumReadTrimmer(14);
    assertEquals(4, bsrt.getTrimPosition(quals, 6));
    bsrt = new BestSumReadTrimmer(13);
    assertEquals(6, bsrt.getTrimPosition(quals, 6));
  }

  public void testGoodBadGoodOk() {
    final byte[] quals = {73, 73, 73, 73, 13, 13, 13, 13, 20, 20, 20, 20};
    BestSumReadTrimmer bsrt = new BestSumReadTrimmer(15);

    assertEquals(quals.length, bsrt.getTrimPosition(quals, quals.length));
    assertEquals(quals.length - 1, bsrt.getTrimPosition(quals, quals.length - 1));
    assertEquals(quals.length - 2, bsrt.getTrimPosition(quals, quals.length - 2));
    assertEquals(quals.length - 3, bsrt.getTrimPosition(quals, quals.length - 3));  //passes because the last value is gt threshold
    assertEquals(4, bsrt.getTrimPosition(quals, 6));
    bsrt = new BestSumReadTrimmer(1);
    assertEquals(quals.length, bsrt.getTrimPosition(quals, quals.length));
    assertEquals(8, bsrt.getTrimPosition(quals, 8));
    bsrt = new BestSumReadTrimmer(14);
    assertEquals(4, bsrt.getTrimPosition(quals, 6));
    bsrt = new BestSumReadTrimmer(13);
    assertEquals(6, bsrt.getTrimPosition(quals, 6));
    bsrt = new BestSumReadTrimmer(23);
    assertEquals(4, bsrt.getTrimPosition(quals, 8));
  }

  public void testZeroLength() {
    final BestSumReadTrimmer bsrt = new BestSumReadTrimmer(65535);
    bsrt.getTrimPosition(new byte[0], 1);
  }

}
