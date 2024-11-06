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

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.alignment.AlignmentResult;

/**
 * Structure holding relevant information for a single read.
 * @param <T> type of hit info subclass.
 */
@TestClass(value = "com.rtg.pairedend.SlidingWindowCollectorTest")
public abstract class AbstractHitInfo<T extends AbstractHitInfo<T>> {
  boolean mFirst;
  boolean mReverseComplement;
  int mReadId;
  int mTemplateStart;
  int mAlignmentScore; // For caching alignment score
  AlignmentResult mAlignment; // For additionally caching a full alignment that was "good enough" to result in a pair being output

  // we keep a linked list of hits, sorted by increasing mTemplateStart.
  private T mNext;
  private T mPrev;

  /**
   * @param first true if hit is from &quot;left&quot; read arm. That is, first in sequencing
   * @param reverseComplement true if read hit on reverse frame
   * @param readId the read ID
   * @param templateStart position of the hit
   */
  void setValues(final boolean first, final boolean reverseComplement, int readId, final int templateStart) {
    mFirst = first;
    mReverseComplement = reverseComplement;
    mReadId = readId;
    mTemplateStart = templateStart;
    mAlignmentScore = -1;
    mAlignment = null;
    mNext = null;
    mPrev = null;
  }

  /**
   * Set the alignment information
   * @param score the alignment score
   * @param alignment the alignment result
   */
  public void setAlignment(int score, AlignmentResult alignment) {
    mAlignmentScore = score;
    mAlignment = alignment;
  }

  /**
   * Insert <code>nextHit</code> after this one.
   *
   * @param nextHit must be non-null and have template start less than or equal to this one.
   */
  void insertHit(T nextHit) {
    link(nextHit);
    mNext = nextHit;
  }

  @SuppressWarnings("unchecked")
  private T thisT() {
    return (T) this;
  }

  /**
   * By separating this part out it prevents the necessity of an unchecked cast
   * to <code>T</code> when assigning to <code>mNext</code> while still allowing
   * <code>mNext</code> to be private.
   */
  private void link(AbstractHitInfo<T> nextHit) {
    assert nextHit.mNext == null;
    assert nextHit.mPrev == null;
    assert nextHit.mTemplateStart >= mTemplateStart;
    nextHit.mNext = mNext;
    if (mNext != null) {
      mNext.setPrev(nextHit.thisT());
    }
    nextHit.mPrev = thisT();
  }

  void setNext(T nextHit) {
    mNext = nextHit;
  }

  void setPrev(T prev) {
    mPrev = prev;
  }

  /**
   * @return the alignment result associated with this hit
   */
  public AlignmentResult alignment() {
    return mAlignment;
  }

  /**
   * @return the alignment score for this hit
   */
  public int score() {
    return mAlignmentScore;
  }

  /**
   * @return true if this hit is from the read that was first in sequencing
   */
  public boolean first() {
    return mFirst;
  }

  /**
   * @return true if this hit is reverse complement
   */
  public boolean reverseComplement() {
    return mReverseComplement;
  }

  /**
   * @return the id of the read for this hit
   */
  public int readId() {
    return mReadId;
  }

  /**
   * @return the template start position for this hit
   */
  public int templateStart() {
    return mTemplateStart;
  }

  T next() {
    return mNext;
  }

  T prev() {
    return mPrev;
  }

  boolean isPair(AbstractHitInfo<T> other) {
    return other != null && mReadId == other.mReadId && mFirst != other.mFirst;
  }
}
