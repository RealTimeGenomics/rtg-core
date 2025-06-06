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
import com.rtg.ngs.NgsParams;
import com.rtg.ngs.NgsParamsBuilder;


/**
 * Tests corresponding class
 */
public class NoIndelsEditDistanceTest extends AbstractUnidirectionalEditDistanceTest {

  @Override
  protected UnidirectionalEditDistance getEditDistanceInstance(int gapOpenPenalty, int gapExtendPenalty, int subsPenalty, int unknownsPenalty) {
    final NgsParams params = new NgsParamsBuilder().gapOpenPenalty(gapOpenPenalty).gapExtendPenalty(gapExtendPenalty).substitutionPenalty(subsPenalty).unknownsPenalty(unknownsPenalty).create();
    return new NoIndelsEditDistance(params);
  }

  public void alignWithoutIndels(String read, String template, String actionsString, int zeroBasedPos, int mismatches, int unknownsPenalty, int subsPenalty) {
    final NgsParams params = new NgsParamsBuilder().gapOpenPenalty(1).gapExtendPenalty(1).substitutionPenalty(subsPenalty).unknownsPenalty(unknownsPenalty).create();
    final NoIndelsEditDistance noied = new NoIndelsEditDistance(params);
    final byte[] encodedRead = DnaUtils.encodeString(read);
    final byte[] encodedTemplate = DnaUtils.encodeString(template);
    final int[] actions = noied.calculateEditDistance(encodedRead, encodedRead.length, encodedTemplate, zeroBasedPos, Integer.MAX_VALUE, 1, true);
    if (mismatches < 0) {
      assertNull(actions);
    } else {
      assertEquals(mismatches, read.length() - ActionsHelper.matchCount(actions));
      assertEquals(mismatches * subsPenalty, ActionsHelper.alignmentScore(actions));
      assertEquals(zeroBasedPos, ActionsHelper.zeroBasedTemplateStart(actions));
      assertEquals(actionsString, ActionsHelper.toString(actions));
    }
  }

  public void alignWithoutIndels(String read, String template, String actionsString, int zeroBasedPos, int mismatches, boolean treatNsAsMismatches) {
    alignWithoutIndels(read, template, actionsString, zeroBasedPos, mismatches, treatNsAsMismatches ? 1 : 0, 1);
  }

  /**
   * To test alignment speed, set repeatCount above to 100,000
   * then run this test.
   */
  public void testAlignWithoutIndelsSpeed36() {
    alignWithoutIndels("atgtacgtacgtacgtacgtacgtacgtacgtacgt",
        "acgtacgtacgtacgtacgtacgtacgtacgtacgt",
        "=X==================================", 0, 1, false);
  }

  /**
   * To test alignment speed on a longer sequence, set repeatCount above
   * to 100,000 then run this test.
   */
  public void testAlignWithoutIndelsSpeed130() {
    final String r32 = "aaaaaaaaccccccccggggggggtttttttt";
    alignWithoutIndels(r32 + r32 + r32 + "ag" + r32,
        r32 + r32 + r32 + "cg" + r32,
        "================================"
            + "================================"
            + "================================"
            + "X="
            + "================================", 0, 1, false);
  }

  public void testAlignWithoutIndels1() {
    alignWithoutIndels("a", "a", "=", 0, 0, false);
  }

  public void testAlignWithoutIndelsMismatchAtStart() {
    alignWithoutIndels("acgt", "tcgt", "X===", 0, 1, false);
  }

  public void testResizeWorkspace() {
    final NgsParams params = new NgsParamsBuilder().gapOpenPenalty(1).gapExtendPenalty(1).substitutionPenalty(1).create();
    final NoIndelsEditDistance noied = new NoIndelsEditDistance(params);

    final byte[] encodedRead = DnaUtils.encodeString("acgt");
    final byte[] encodedTemplate = DnaUtils.encodeString("tcgttcgttgtgtg");
    final int[] actions = noied.calculateEditDistance(encodedRead, encodedRead.length, encodedTemplate, 0, Integer.MAX_VALUE, 1, true);
    assertNotNull(actions);
    assertEquals(1, encodedRead.length - ActionsHelper.matchCount(actions));
    assertEquals(1, ActionsHelper.alignmentScore(actions));
    assertEquals(0, ActionsHelper.zeroBasedTemplateStart(actions));

    final byte[] encodedReadLonger = DnaUtils.encodeString("gtacgtt");
    final int[] actions2 = noied.calculateEditDistance(encodedReadLonger, encodedReadLonger.length, encodedTemplate, 2, Integer.MAX_VALUE, 1, true);
    assertNotNull(actions2);
    assertEquals(1, encodedReadLonger.length - ActionsHelper.matchCount(actions2));
    assertEquals(1, ActionsHelper.alignmentScore(actions2));
    assertEquals(2, ActionsHelper.zeroBasedTemplateStart(actions2));
  }

  public void testAlignWithoutIndelsMismatchAtEnd() {
    alignWithoutIndels("acgt", "cacgacc", "===X", 1, 1, false);
  }

  public void testAlignWithoutIndelsTooManyMismatches() {
    //default penalties of 1,1,1 so need 3 mismatches to stop.
    alignWithoutIndels("acgttcta", "tcgatcga", "", 0, -1, false);
  }

  public void testAlignWithoutIndelsTemplateMuchTooShort() {
    alignWithoutIndels("act", "act", "", 2, -1, false);
  }

  public void testAlignWithoutIndelsTemplateTooShort() {
    alignWithoutIndels("act", "act", "", 1, -1, false);
  }

  public void testAlignWithoutIndelsAllNs() {
    alignWithoutIndels("nnnn", "nnnn", "====", 0, -1, false);
  }

  public void testAlignWithoutIndelsMostlyNs() {
    alignWithoutIndels("ncnn", "nntn", "====", 0, -1, false);
  }

  public void testAlignWithoutIndelsTemplate1N() {
    alignWithoutIndels("tttt", "cnttt", "====", 1, -1, false);
  }

  public void testAlignWithoutIndelsEndMismatch() {
    alignWithoutIndels("tggac", "tggag", "====X", 0, 1, false);
  }

  public void testAlignWithoutIndelsFirstMismatch() {
    alignWithoutIndels("aggac", "tggac", "X====", 0, 1, false);
  }

  public void testAlignWithoutIndelsEndTwoMismatches() {
    alignWithoutIndels("tggac", "tggtg", "===XX", 0, 2, false);
  }

  public void testAlignWithoutIndelsStartTwoMismatches() {
    alignWithoutIndels("acgtg", "tggtg", "XX===", 0, 2, false);
  }
}
