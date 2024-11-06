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

