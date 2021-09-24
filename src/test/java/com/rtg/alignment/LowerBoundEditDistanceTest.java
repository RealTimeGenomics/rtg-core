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
