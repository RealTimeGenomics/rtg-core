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
package com.rtg.util;

import java.util.Arrays;

import junit.framework.TestCase;

/**
 * Tests corresponding class
 */
public class PortableRandomTest extends TestCase {

  public void testBasics() {

    final PortableRandom pr = new PortableRandom(3);

    assertEquals(-1155099828, pr.nextInt());
    assertEquals(-1879439976, pr.nextInt());
    assertEquals(304908421, pr.nextInt());

    assertEquals(16, pr.nextInt(45));
    assertEquals(8, pr.nextInt(10));
    assertEquals(0, pr.nextInt(1));

    assertEquals(0.768156984078079, pr.nextDouble());
    assertEquals(0.22733466107144407, pr.nextDouble());
    assertEquals(0.6603196166875382, pr.nextDouble());

    assertEquals(-3566243377172859107L, pr.nextLong());
    assertEquals(550039120288444364L, pr.nextLong());
    assertEquals(-3483296404361882349L, pr.nextLong());

    assertTrue(pr.nextBoolean());
    assertTrue(pr.nextBoolean());
    assertFalse(pr.nextBoolean());

    assertEquals(1.6791616586691243, pr.nextGaussian());
    assertEquals(0.6070711938870896, pr.nextGaussian());
    assertEquals(-1.3490545500116493, pr.nextGaussian());

    final byte[] b1 = new byte[5];

    pr.nextBytes(b1);

    assertTrue(Arrays.equals(new byte[] {(byte) 0x51, (byte) 0xD5, (byte) 0xA5, (byte) 0x76, (byte) 0x00}, b1));
    pr.nextBytes(b1);
    assertTrue(Arrays.equals(new byte[] {(byte) 0xF2, (byte) 0xFF, (byte) 0x64, (byte) 0xCF, (byte) 0xF3}, b1));
    pr.nextBytes(b1);
    assertTrue(Arrays.equals(new byte[] {(byte) 0x8C, (byte) 0xDB, (byte) 0xE7, (byte) 0x67, (byte) 0x12}, b1));
  }
}
