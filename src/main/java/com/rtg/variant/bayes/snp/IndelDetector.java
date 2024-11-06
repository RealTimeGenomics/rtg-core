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

package com.rtg.variant.bayes.snp;

import com.rtg.launcher.globals.CoreGlobalFlags;
import com.rtg.launcher.globals.GlobalFlags;
import com.rtg.variant.bayes.EvidenceAcceptor;
import com.rtg.variant.bayes.EvidenceInterface;
import com.rtg.variant.bayes.Model;

/**
 * This class essentially just counts the number of indels seen at a particular point, in order
 * to contribute toward triggering the complex caller.
 */
public class IndelDetector implements EvidenceAcceptor {

  /** The minimum number of indels to occur at a position before invoking the complex caller */
  public static final int TRIVIAL_INDEL_COUNT = 2; //Integer.parseInt(System.getProperty("rtg.modelindel.trivial-indels", "2"));

  static final int MAX_INDEL_EXTENSION = GlobalFlags.getIntegerValue(CoreGlobalFlags.COMPLEX_REGION_SIMPLE_REPEAT_LIMIT);

  private int mNonTrivialInsertCount;
  private int mNonTrivialDeletionCount;
  private int mSoftClipLeftCount;
  private int mSoftClipRightCount;
  private int mNonIndelCount;

  private int mIndelLength;
  private int mSoftClipLength;

  @Override
  public void increment(EvidenceInterface evidence) {
    if (evidence == null) {
      ++mNonIndelCount;
    } else {
      assert evidence instanceof EvidenceIndel;
      if (evidence.mapError() < Model.AMBIGUITY_THRESHOLD) {
        final int opLength = Math.min(((EvidenceIndel) evidence).operationLength(), MAX_INDEL_EXTENSION);
        switch (evidence.read()) {
          case EvidenceIndel.INSERT:
            ++mNonTrivialInsertCount;
            mIndelLength += opLength;
            break;
          case EvidenceIndel.DELETE:
            ++mNonTrivialDeletionCount;
            mIndelLength += opLength;
            break;
          case EvidenceIndel.SOFT_CLIP_LEFT:
            ++mSoftClipLeftCount;
            mSoftClipLength = Math.max(opLength, mSoftClipLength);
            break;
          case EvidenceIndel.SOFT_CLIP_RIGHT:
            ++mSoftClipRightCount;
            mSoftClipLength = Math.max(opLength, mSoftClipLength);
            break;
          default:
            throw new RuntimeException("Invalid indel code: " + evidence.read());
        }
      }
    }
  }

  @Override
  public void freeze() {
    // Nothing to do here
  }

  /**
   * @return typical length of an insert or delete at this position as specified by cigar
   */
  public int indelLength() {
    final int tot = mNonTrivialDeletionCount + mNonTrivialInsertCount;
    return tot == 0 ? 0 : mIndelLength / tot;
  }

  /**
   *
   * @return maximum length of a soft clip as specified by cigar
   */
  public int maxSoftClipLength() {
    return mSoftClipLength;
  }

  /**
   * Get non-trivial insert count.
   * @return Returns the non-trivial insert count.
   */
  int nonTrivialInsertCount() {
    return mNonTrivialInsertCount;
  }

  /**
   * Get non-trivial deletion count.
   * @return Returns the non-trivial deletion count.
   */
  int nonTrivialDeletionCount() {
    return mNonTrivialDeletionCount;
  }

  int softClipLeftCount() {
    return mSoftClipLeftCount;
  }

  int softClipRightCount() {
    return mSoftClipRightCount;
  }

  int minIndelCount(double fraction) {
    return Math.max(TRIVIAL_INDEL_COUNT, (int) (mNonIndelCount * fraction));
  }

  int totalCount() {
    return mNonIndelCount + mNonTrivialInsertCount + mNonTrivialDeletionCount + mSoftClipLeftCount + mSoftClipRightCount;
  }
}
