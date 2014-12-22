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

  /**
   * @param first true if hit is from &quot;left&quot; read arm. That is, first in sequencing
   */
  void setValues(final boolean first, final boolean reverseComplement, int readId, final int templateStart) {
    mFirst = first;
    mReverseComplement = reverseComplement;
    mReadId = readId;
    mTemplateStart = templateStart;
    mAlignmentScore = -1;
    mAlignment = null;
    mNext = null;
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

  /**
   * By separating this part out it prevents the necessity of an unchecked cast
   * to <code>T</code> when assigning to <code>mNext</code> while still allowing
   * <code>mNext</code> to be private.
   */
  private void link(AbstractHitInfo<T> nextHit) {
    assert nextHit.mNext == null;
    assert nextHit.mTemplateStart >= mTemplateStart;
    nextHit.mNext = mNext;
  }

  void setNext(T nextHit) {
    mNext = nextHit;
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

  boolean isPair(AbstractHitInfo<T> other) {
    return other != null && mReadId == other.mReadId && mFirst != other.mFirst;
  }
}
