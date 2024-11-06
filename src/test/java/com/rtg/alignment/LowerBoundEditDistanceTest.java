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
package com.rtg.alignment;

import com.rtg.mode.DnaUtils;

/**
 */
public class LowerBoundEditDistanceTest extends AbstractUnidirectionalEditDistanceTest {

  @Override
  protected UnidirectionalEditDistance getEditDistanceInstance(int gapOpenPenalty, int gapExtendPenalty, int subsPenalty, int unknownsPenalty) {
    return new LowerBoundEditDistance(4, 1, unknownsPenalty);
  }

  /**
   * Simple test
   */
  public void testEmpty() {
    final byte[] read = DnaUtils.encodeString("aaaaaaaaa");
    final byte[] temp = DnaUtils.encodeString("ttttttttt");
    final UnidirectionalEditDistance lo = new LowerBoundEditDistance(4, 1, 1);
    final int[] res = lo.calculateEditDistance(read, read.length, temp, 1, 0, 7, true);
    assertNotNull(res);
   }


  /**
   * Simple test
   */
  public void testEmpty2() {
    final byte[] read = DnaUtils.encodeString("aaaaaaaaa");
    final byte[] temp = DnaUtils.encodeString("ttttttttt");
    final UnidirectionalEditDistance lo = new LowerBoundEditDistance(4, 1, 1);
    final int[] res = lo.calculateEditDistance(read, read.length, temp, 0, 2, 7, true);
    assertNull(res); // null means it's under the threshold
  }

  /**
   * Simple test
   */
  public void testEmpty3() {
    final byte[] read = DnaUtils.encodeString("aaaaaaaaa");
    final byte[] temp = DnaUtils.encodeString("ttttttttt");
    final UnidirectionalEditDistance lo = new LowerBoundEditDistance(4, 1, 1);
    final int[] res = lo.calculateEditDistance(read, read.length, temp, 0, 0, 7, true);
    assertEquals(Integer.MAX_VALUE, res[ActionsHelper.ALIGNMENT_SCORE_INDEX]);
    assertEquals(ActionsHelper.ACTIONS_START_INDEX, res.length);
  }

  /**
   * Simple test
   */
  public void testEmpty4() {
    final byte[] read = DnaUtils.encodeString("aaaaaaaaa");
    final byte[] temp = DnaUtils.encodeString("ttttttttt");
    final UnidirectionalEditDistance lo = new LowerBoundEditDistance(4, 1, 1);
    final int[] res = lo.calculateEditDistance(read, read.length, temp, 0, Integer.MAX_VALUE, 7, true);
    assertNull(res); // null means under the threshold
  }

  /**
   * Simple test
   */
  public void test5() {
    final byte[] read = DnaUtils.encodeString("GGAGGCCCTCTTGTCTTGCA");
    final byte[] temp = DnaUtils.encodeString("TGCAAGACAAGAGGGCCTCC");
    final UnidirectionalEditDistance lo = new LowerBoundEditDistance(4, 1, 1);
    final int[] res = lo.calculateEditDistance(read, read.length, temp, 0, 1, 7, true);
    assertEquals(Integer.MAX_VALUE, res[ActionsHelper.ALIGNMENT_SCORE_INDEX]);
  }

  /**
   * Simple test
   */
  public void test5r() {
    final byte[] read = DnaUtils.encodeString("TGCAAGACAAGAGGGCCTCC");
    final byte[] temp = DnaUtils.encodeString("GGAGGCCCTCTTGTCTTGCA");
    final UnidirectionalEditDistance lo = new LowerBoundEditDistance(4, 1, 1);
    final int[] res = lo.calculateEditDistance(read, read.length, temp, 0, 2, 7, true);
    assertEquals(Integer.MAX_VALUE, res[ActionsHelper.ALIGNMENT_SCORE_INDEX]);
  }


  /**
   * Simple test
   */
  public void test6() {
    final byte[] read = DnaUtils.encodeString("GGAGGCCCTCTTGTCTTGCANN");
    final byte[] temp = DnaUtils.encodeString("TGCAAGACAAGAGGGCCTCC");
    final UnidirectionalEditDistance lo = new LowerBoundEditDistance(4, 1, 1);
    final int[] res = lo.calculateEditDistance(read, read.length, temp, 7, 4, 7, true);
    assertEquals(Integer.MAX_VALUE, res[ActionsHelper.ALIGNMENT_SCORE_INDEX]);
  }


}
