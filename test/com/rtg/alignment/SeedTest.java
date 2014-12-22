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

package com.rtg.alignment;

import com.rtg.util.Utils;

import junit.framework.TestCase;

/**
 */
public class SeedTest extends TestCase {

  public void test() {
    final Seed sd = new Seed(4);
    sd.integrity();
    assertEquals("size=4 mask=11111111:00000000:00000000:11111111", sd.toString());
    assertEquals("11111111:00000000:00000000:11111111", Utils.toBitsSep(sd.init()));
    assertEquals(4, sd.size());

    assertEquals(0, sd.next(0, (byte) 1));
    assertEquals(1, sd.next(0, (byte) 2));
    assertEquals(2, sd.next(0, (byte) 3));
    assertEquals(3, sd.next(0, (byte) 4));
    assertEquals("11111111:00000000:00000000:11111111", Utils.toBitsSep(sd.next(0, (byte) 0)));

    int seed = 0;
    for (byte b = 1; b <= 4; b++) {
      seed = sd.next(seed, b);
      assertTrue(sd.isValid(seed));
    }
    assertEquals("00000000:00000000:00000000:00011011", Utils.toBitsSep(seed));

    seed = sd.next(seed, (byte) 0);
    assertFalse(sd.isValid(seed));
    assertEquals("11111111:00000000:00000000:11111111", Utils.toBitsSep(seed));

    seed = sd.next(seed, (byte) 3);
    assertFalse(sd.isValid(seed));
    assertEquals("11111100:00000000:00000000:11111110", Utils.toBitsSep(seed));

    seed = sd.next(seed, (byte) 3);
    assertFalse(sd.isValid(seed));
    assertEquals("11110000:00000000:00000000:11111010", Utils.toBitsSep(seed));

    seed = sd.next(seed, (byte) 3);
    assertFalse(sd.isValid(seed));
    assertEquals("11000000:00000000:00000000:11101010", Utils.toBitsSep(seed));

    seed = sd.next(seed, (byte) 3);
    assertTrue(sd.isValid(seed));
    assertEquals("00000000:00000000:00000000:10101010", Utils.toBitsSep(seed));
  }

  public void testBad() {
    try {
      new Seed(8);
      fail();
    } catch (final Exception e) {
      //expected
    }
  }
}
