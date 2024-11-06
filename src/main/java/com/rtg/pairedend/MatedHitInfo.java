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
package com.rtg.pairedend;


import com.rtg.alignment.AlignmentResult;

/**
 * A simple structure to store information about a pair of hits.
 * Left/right refers to placement on the reference.
 * First/second refers to whether the read arm was first or second in sequencing.
 *
 */
public class MatedHitInfo {

  int mReadId;
  boolean mFirstRight;
  boolean mReverseComplementRight;
  int mTemplateStartRight;
  boolean mReverseComplementLeft;
  int mTemplateStartLeft;

  // Populated and stored here so that alignments for a mating need only be calculated once
  AlignmentResult mAlignmentLeft = null;
  int mAlignmentScoreLeft = -1;
  AlignmentResult mAlignmentRight = null;
  int mAlignmentScoreRight = -1;

  // Used to indicate that the left alignment is so poor that there is no point attempting further matings with this hit
  boolean mLeftIsVeryPoor;

  /**
   * Helper method to set many values at once
   * @param readId the read id
   * @param firstRight true if the read that was first in sequencing is rightmost on the reference
   * @param reverseComplementLeft left read is reverse complement
   * @param templateStartLeft template start of the leftmost read
   * @param reverseComplementRight if the right is reverse complement
   * @param templateStartRight template start of the rightmost read
   */
  public void setValues(int readId, boolean firstRight, boolean reverseComplementLeft, int templateStartLeft, boolean reverseComplementRight, int templateStartRight) {
    mReadId = readId;
    mFirstRight = firstRight;
    mReverseComplementRight = reverseComplementRight;
    mTemplateStartRight = templateStartRight;
    mReverseComplementLeft = reverseComplementLeft;
    mTemplateStartLeft = templateStartLeft;
    mLeftIsVeryPoor = false;
    setAlignments(null, null);
  }

  /**
   * Set data about this hit
   * @param left the alignment result of the leftmost alignment
   * @param right the alignment result of the rightmost alignment
   */
  public void setAlignments(AlignmentResult left, AlignmentResult right) {
    mAlignmentLeft = left;
    mAlignmentScoreLeft = mAlignmentLeft == null ? -1 : left.getScore();
    mAlignmentRight = right;
    mAlignmentScoreRight = mAlignmentRight == null ? -1 : right.getScore();
  }

  public int getAlignmentScoreLeft() {
    return mAlignmentScoreLeft;
  }

  public void setAlignmentScoreLeft(int alignmentScoreLeft) {
    mAlignmentScoreLeft = alignmentScoreLeft;
  }

  public int getAlignmentScoreRight() {
    return mAlignmentScoreRight;
  }

  public void setAlignmentScoreRight(int alignmentScoreRight) {
    mAlignmentScoreRight = alignmentScoreRight;
  }

  public void setLeftVeryPoor(boolean leftIsVeryPoor) {
    mLeftIsVeryPoor = leftIsVeryPoor;
  }

  public int getReadId() {
    return mReadId;
  }

  public AlignmentResult getAlignmentLeft() {
    return mAlignmentLeft;
  }

  public AlignmentResult getAlignmentRight() {
    return mAlignmentRight;
  }

  public boolean isFirstRight() {
    return mFirstRight;
  }

  public int getTemplateStartLeft() {
    return mTemplateStartLeft;
  }

  public boolean isReverseComplementRight() {
    return mReverseComplementRight;
  }

  public int getTemplateStartRight() {
    return mTemplateStartRight;
  }

  public boolean isReverseComplementLeft() {
    return mReverseComplementLeft;
  }
}

