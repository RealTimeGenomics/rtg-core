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
package com.rtg.index.hash;

import com.rtg.mode.DNA;

/**
 */
public class ExactHashFunctionTest extends AbstractHashFunctionTest {

  @Override
  protected HashFunction getHashFunction(final int windowSize, final int bits) {
    return new ExactHashFunction(windowSize, bits);
  }

  /**
   * Test method for {@link com.rtg.index.hash.ExactHashFunction}.
   */
  public final void test1() {
    final char[] codes = DNA.valueChars();
    final ExactHashFunction hf = (ExactHashFunction) getHashFunction(3, 2);
    assertFalse(hf.isValid());
    hf.integrity();
    assertEquals("AAA", hf.hashToSeq(0, codes));

    assertEquals(0, hf.hashStep((byte) 0));
    assertFalse(hf.isValid());

    assertEquals(1, hf.hashStep((byte) 1));
    assertFalse(hf.isValid());
    assertEquals("AAC", hf.hashToSeq(1, codes));

    assertEquals(5, hf.hashStep((byte) 1));
    assertTrue(hf.isValid());
    assertEquals("ACC", hf.hashToSeq(5, codes));

    assertEquals(21, hf.hashStep((byte) 1));
    assertTrue(hf.isValid());
    assertEquals("CCC", hf.hashToSeq(21, codes));

    assertEquals(21, hf.hashStep((byte) 1));
    assertTrue(hf.isValid());

    assertEquals(22, hf.hashStep((byte) 2));
    assertTrue(hf.isValid());
    assertEquals("CCG", hf.hashToSeq(22, codes));

    assertEquals(27, hf.hashStep((byte) 3));
    assertTrue(hf.isValid());
    assertEquals("CGT", hf.hashToSeq(27, codes));
    assertEquals("TTT", hf.hashToSeq(63, codes));

    hf.reset();
    assertFalse(hf.isValid());
    hf.integrity();
    assertEquals(0, hf.hashStep((byte) 0));
    assertFalse(hf.isValid());
    assertEquals(1, hf.hashStep((byte) 1));
    assertFalse(hf.isValid());
  }

  public void testDualMode() {
    final char[] codes = DNA.valueChars();
    final ExactHashFunction hf = new ExactHashFunction(3, 2, true);
    assertFalse(hf.isValid());
    hf.integrity();
    //0 -> A, 1 -> C, 2 -> G, 3 ->T
    //we are testing reverse complement
    assertEquals(48, hf.reverseHashStep((byte) 0));
    assertEquals("TAA", hf.hashToSeq(48, codes));
    assertEquals(44, hf.reverseHashStep((byte) 1));
    assertEquals("GTA", hf.hashToSeq(44, codes));
    assertEquals(43, hf.reverseHashStep((byte) 1));
    assertEquals("GGT", hf.hashToSeq(43, codes));
    assertEquals(42, hf.reverseHashStep((byte) 1));
    assertEquals("GGG", hf.hashToSeq(42, codes));
    assertEquals(42, hf.reverseHashStep((byte) 1));
    assertEquals("GGG", hf.hashToSeq(42, codes));
    assertEquals(26, hf.reverseHashStep((byte) 2));
    assertEquals("CGG", hf.hashToSeq(26, codes));
    assertEquals(6, hf.reverseHashStep((byte) 3));
    assertEquals("ACG", hf.hashToSeq(6, codes));
  }

  /**
   * Test method for {@link com.rtg.util.integrity.IntegralAbstract#toString()}.
   */
  public final void testToString() {
    final ExactHashFunction hf = (ExactHashFunction) getHashFunction(3, 2);
    hf.integrity();
    assertEquals("SimpleHashFunction window size=3 bits=2 mask=63 hash=0", hf.toString());
    assertEquals(3, hf.hashStep((byte) 3));
    assertEquals("SimpleHashFunction window size=3 bits=2 mask=63 hash=3", hf.toString());
  }

}

