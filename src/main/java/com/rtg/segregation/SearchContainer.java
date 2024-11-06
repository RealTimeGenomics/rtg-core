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

package com.rtg.segregation;

import com.rtg.util.Pair;
import com.rtg.util.Utils;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;

/**
 * Holds each entry of a search path.
 */
class SearchContainer extends IntegralAbstract implements Comparable<SearchContainer> {

  /** Previous container in the search chain. null at start. */
  private final SearchContainer mPrevLink;
  /** Cumulative score for chain. Smaller is better. */
  private final double mScore;
  /** Last block that was phase compatible. null if this block is phase compatible. */
  private final SearchContainer mLastOK;
  /** The block corresponding to this container (that is, at this position). */
  private final SegregationBlock mBlock;
  /** The phasing pattern for this block and for previous blocks that are compatible (may be more precise than them). */
  private final PatternArray mPattern;
  /** The type of block - whether it is phase compatible or not etc. */
  private final SearchType mType;
  /** If this is the block immediately following a cross over then a non-null cross over description, otherwise null. */
  private final CrossOver mXO;
  /** An identifier unique within any one sequence of a search. It is generated externally and is used to provide a unique compare to and allow debugging. */
  private final long mId;


  /**
   * @param prevLink previous link in chain - can be null at start of chain.
   * @param score smaller is better.
   * @param lastOK last container that is phase compatible.
   * @param block the current block.
   * @param pattern the pattern for this block
   * @param type of the container.
   * @param xo cross over - null unless search type is <code>XO</code>.
   * @param id unique identifier
   */
  SearchContainer(SearchContainer prevLink, double score, SearchContainer lastOK, SegregationBlock block, PatternArray pattern, SearchType type, CrossOver xo, long id) {
    mPrevLink = prevLink;
    mScore = score;
    mLastOK = lastOK;
    mBlock = block;
    mPattern = pattern;
    mType = type;
    mXO = xo;
    mId = id;
    assert integrity();
  }

  /**
   * @return the cumulative score of this search path.
   */
  double score() {
    return mScore;
  }

  /** @return the pattern of this block and any previous phase compatible containers. May be more precise than individual blocks. */
  PatternArray pattern() {
    return mPattern;
  }

  /**
   * @return a unique key used when working out if two paths have converged back to the same state.
   */
  Pair<PatternArray, Boolean> key() {
    return new Pair<>(mPattern, mType == SearchType.Error);
  }

  /**
   * @return the last container that was phase compatible. May be this one.
   */
  SearchContainer goodContainer() {
    if (mType == SearchType.Error) {
      return mLastOK;
    }
    return this;
  }

  /**
   * @return the current block.
   */
  SegregationBlock block() {
    return mBlock;
  }

  /**
   * @return previous container in chain. null at start.
   */
  SearchContainer previous() {
    return mPrevLink;
  }

  /**
   * @return the type of this container.
   */
  SearchType type() {
    return mType;
  }

  /**
   * @return if this is a container that immediately follows a crossover then the cross over descriptor, otherwise null.
   */
  CrossOver xo() {
    return mXO;
  }

  @Override
  public int compareTo(SearchContainer that) {
    final int c = Double.compare(mScore, that.mScore);
    if (c != 0) {
      return c;
    }
    return (int) (this.mId - that.mId);
  }

  @Override
  public int hashCode() {
    return (int) (mId ^ (mId >> 32));
  }

  @Override
  public boolean equals(Object obj) {
    if (obj == this) {
      return true;
    }
    if (obj == null) {
      return false;
    }
    final SearchContainer that = (SearchContainer) obj;
    return this.mId == that.mId;
  }

  @Override
  public void toString(StringBuilder sb) {
    sb.append(mId).append(".").append(mType);
    sb.append("<").append(mPrevLink == null ? "null" : mPrevLink.mId);
    sb.append(":").append(Utils.realFormat(mScore, 1)).append(" fa: ").append(mPattern.faString()).append(" mo: ").append(mPattern.moString());
  }

  @Override
  public final boolean integrity() {
    Exam.assertTrue(Double.isFinite(mScore));
    Exam.assertNotNull(mBlock);
    Exam.assertNotNull(mPattern);
    Exam.assertNotNull(mType);
    if (mLastOK != null) {
      Exam.assertNotNull(mPrevLink);
    }
    if (mPrevLink == null) {
      Exam.assertTrue(SearchType.New == mType || SearchType.Error == mType);
    }
    switch (mType) {
      case OK:
        Exam.assertTrue(mXO == null);
        Exam.assertTrue(mLastOK == null);
        Exam.assertTrue(mPattern.compatible(mBlock.patterns()));
        break;
      case Error:
        Exam.assertTrue(mXO == null);
        break;
      case XO:
        Exam.assertNotNull(mXO);
        Exam.assertTrue(mLastOK == null);
        //Exam.assertTrue(mPattern.strictEquals(mBlock.patterns()));
        break;
      case New:
        Exam.assertTrue(mXO == null);
        Exam.assertTrue(mLastOK == null);
        Exam.assertTrue(mPattern.strictEquals(mBlock.patterns()));
        break;
      default:
        throw new RuntimeException();
    }
    return true;
  }

}
