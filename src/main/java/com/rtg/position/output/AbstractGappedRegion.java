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
package com.rtg.position.output;

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.launcher.BuildParams;
import com.rtg.mode.Frame;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;

/**
 * Holds a HashingRegion which consists of one or more exact match word
 * hits separated by gaps.
 * @param <G> type of subclass (needed to get some of the parameters and return results correct).
 */
@TestClass("com.rtg.position.output.GappedRegionTest")
public abstract class AbstractGappedRegion<G extends AbstractGappedRegion<G>> extends IntegralAbstract {

  enum State {
    EMPTY, INIT, WRITTEN
  }

  protected static final int NULL = -1;

  private long mBucket = NULL;

  final int mId;

  protected final int mWordSize;

  protected final int mStepSize;

  protected final GapScorer mDist;

  protected final Frame[] mSubjectFrames;

  State mState = State.EMPTY;

  protected int mBuildSeqId = NULL;

  int mBuildStart = NULL;

  int mBuildEnd = NULL;

  int mQueryStart = NULL;

  int mQueryEnd = NULL;

  boolean mReverseFrame;

  private GapBuckets<G> mBuckets = null;

  private G mNext = null;

  /**
   * @param id an identifier unique within a particular <code>GappedOutput</code> object.
   * @param params parameters for the query.
   * @param distribution probability distribution for gap probabilities.
   */
  AbstractGappedRegion(final int id, final BuildParams params, final GapScorer distribution) {
    mId = id;
    mStepSize = params.stepSize();
    mWordSize = params.windowSize();
    mDist = distribution;
    mSubjectFrames = params.sequences().mode().allFrames();
  }

  /**
   * Check if initialized.
   * @return true iff initialized.
   */
  public boolean isInitialized() {
    return mState == State.INIT;
  }

  /**
   * Check if initialized or written (that is structure still valid).
   * @return true iff initialized.
   */
  public boolean isValid() {
    return mState == State.INIT || mState == State.WRITTEN;
  }

  /**
   * Initialize the region.
   * This instead of a constructor so that the set of regions can be all pre-constructed and
   * then we just initialize and reset them.
   * This assumes the initial region is of one word length.
   * @param seqId internal sequence identifier.
   * @param buildStart position of start of build sequence (&gt;= 0).
   * @param queryStart position of start of query sequence (&gt;= 0).
   * @param gb the buckets data structure we interact with via bucket().
   * @param next region in chain from a bucket.
   * @param reverseFrame we switch the start and end of the build when running in reverse frame to allow merging of two regions to
   *        result in a region that remains in the same bucket.
   */
  void initialize(final int seqId, final int buildStart, final int queryStart, final GapBuckets<G> gb, final G next, final boolean reverseFrame) {
    assert mState == State.EMPTY;
    //assert buildStart % mStepSize == 0;
    assert next != null;
    mBuildSeqId = seqId;
    mReverseFrame = reverseFrame;
    //if (reverseFrame) {
    //  mBuildStart = buildStart + mWordSize - 1;
    //  mBuildEnd = buildStart;
    //  mQueryStart = queryStart;
    //  mQueryEnd = queryStart + mWordSize - 1;
      //mQueryStart = queryStart + mWordSize - 1;
      //mQueryEnd = queryStart;
    //} else {
      mBuildStart = buildStart;
      mBuildEnd = buildStart + mWordSize - 1;
      mQueryStart = queryStart;
      mQueryEnd = queryStart + mWordSize - 1;
    //}
    mBuckets = gb;
    mBucket = mBuckets.bucket(mBuildSeqId, mBuildEnd, mQueryEnd);
    mNext = next;
    mState = State.INIT;
    //System.err.println("initialize " + this);
  }

  /**
   * Clear it out so that can later be reused via initialize.
   */
  void reset() {
    assert mState == State.WRITTEN;
    mState = State.EMPTY;
    mBuildSeqId = NULL;
    mBuildStart = NULL;
    mBuildEnd = NULL;
    mQueryStart = NULL;
    mQueryEnd = NULL;
    mBuckets = null;
    mBucket = NULL;
    mNext = null;
  }

  /**
   * Return the probability of the gap between this and that regions.
   * that will usually be the "earlier" (otherwise 0.0 is returned).
   * @param that the other region being compared against.
   * @return the probability.
   */
  double gapScore(final G that) {
    return gapScore(mDist, that);
  }

  /**
   * Return the score of the gap between this and that regions.
   * that will usually be the "earlier" (otherwise NEGATIVE_INFINITY is returned).
   * @param dist the probability distribution for gaps that has been pre-computed.
   * @param that the other region being compared against.
   * @return the score.
   */
  double gapScore(final GapScorer dist, final G that) {
    //System.err.println("gapProbability this=" + this + " that=" + that);
    assert isInitialized();
    if (this == that) {
      return Double.NEGATIVE_INFINITY;
    }
    assert this.mQueryEnd >= that.mQueryEnd : this.mQueryEnd + ":" + that.mQueryEnd;
    assert this != that;
    assert this.mId != that.mId;
    assert this.mStepSize == that.mStepSize && this.mWordSize == that.mWordSize;
    assert this.isInitialized() && that.isInitialized() : "this=" + this + " that=" + that;
    if (this.mBuildSeqId != that.mBuildSeqId) {
      //this allows buckets to be sloppy and contain different sequences
      return Double.NEGATIVE_INFINITY;
    }
    dist.setInternalBuildId(mBuildSeqId);
    //TODO deal with reverse frame
    //if (mReverseFrame) {
    //  return dist.score(this.mBuildStart + 1, that.mBuildEnd, this.mQueryStart + 1, that.mQueryEnd);
    //} else {
      return dist.score(that.mBuildEnd + 1, this.mBuildStart, that.mQueryEnd + 1, this.mQueryStart);
    //}
  }

  /**
   * Merge that region into this region.
   * Reset that for later re-use.
   * that will usually be the "earlier".
   * @param that the region to be merged.
   */
  void merge(final G that) {
    assert this.mBuildSeqId == that.mBuildSeqId;
    //assert this.mBuildStart > that.mBuildStart;
    //assert this.mQueryStart > that.mQueryStart;
    assert this.mBuildEnd > that.mBuildEnd : this.mBuildEnd + ":" + that.mBuildEnd;
    assert this != that;
    assert this.mId != that.mId;
    assert this.mStepSize == that.mStepSize && this.mWordSize == that.mWordSize;
    assert this.isInitialized() && that.isInitialized();
    //TODO deal with reverse frame
    if (this.mBuildStart <= that.mBuildStart) {
      assert (this.mBuildEnd - this.mQueryEnd) == (that.mBuildEnd - that.mQueryEnd) : "this buildEnd: " + this.mBuildEnd + " this queryEnd: " + this.mQueryEnd + " that buildEnd: " + that.mBuildEnd + " that queryEnd: " + that.mQueryEnd + " rev?" + mReverseFrame + " this: " + this + " that: " + that;
      //do nothing the start is already long enough
      return;
    }
    mBuildStart = that.mBuildStart;
    mQueryStart = that.mQueryStart;
  }

  /**
   * Get the bucket.
   * This is computed using the end points of the build and query sequences.
   * @return the bucket.
   */
  long bucket() {
    assert isInitialized();
    return mBucket;
  }

  /**
   * Get the first bucket for the build sequence.
   * This is computed using the end points of the build and query sequences.
   * @return the bucket.
   */
  long firstBucket() {
    assert isInitialized();
    return mBuckets.firstBucket(mBuildSeqId);
  }

  /**
   * Get the last bucket for the build sequence.
   * This is computed using the end points of the build and query sequences.
   * @return the bucket.
   */
  long lastBucket() {
    assert isInitialized();
    return mBuckets.lastBucket(mBuildSeqId);
  }

  /**
   * Get the bucket after applying an offset to the delta.
   * This is computed using the end points of the build and query sequences.
   * @param delta the delta
   * @return the bucket.
   */
  long bucket(final int delta) {
    assert isInitialized();
    return mBuckets.bucket(mBuildSeqId, mBuildEnd + delta, mQueryEnd);
  }

  /**
   * Get a bucket end point after applying an offset to the delta.
   * This is computed using the end points of the build and query sequences.
   * It differs from bucket in that it doesnt wrap at the ends so the
   * returned value may lie outside the first to last range.
   * @param delta the delta
   * @return the bucket.
   */
  long bucketEx(final int delta) {
    assert isInitialized();
    return mBuckets.bucketEx(mBuildSeqId, mBuildEnd + delta, mQueryEnd);
  }

  /**
   * Set the value of next.
   * @param next region to set next to.
   */
  void setNext(final G next) {
    assert isInitialized();
    assert next != null;
    mNext = next;
  }

  /**
   * Get the next region.
   * @return the next region.
   */
  G next() {
    assert isValid();
    return mNext;
  }

  /**
   * Get  the build sequence id.
   * @return the build sequence id.
   */
  int sequenceId() {
    assert isInitialized();
    return mBuildSeqId;
  }

  /**
   * Get the unique id of the region.
   * @return the unique id of the region.
   */
  int id() {
    return mId;
  }

  /**
   * Get  the end position for the query.
   * Used in a number of places for testing.
   * @return the end position. (&gt;= word size &lt; query sequence length)
   */
  int queryEnd() {
    assert isInitialized();
    return mQueryEnd;
  }

  /**
   * Get total length of region (in build co-ordinates may differ from
   * read side in presence of indels).
   * @return total length of region.
   */
  int buildLength() {
    assert isInitialized();
    final int length = mBuildEnd - mBuildStart + 1;
    assert (length - mWordSize) % mStepSize == 0;
    return length;
  }

  protected void write() {
    assert isInitialized() : mState;
    mState = State.WRITTEN;
  }

  SurrogateRegion surrogateOutput(final SurrogateRegion surrogate, final int queryId, final Frame queryFrame, final int queryLength, final int queryEffectiveLength) {
    assert isInitialized() : mState;
    final SurrogateRegion res;
    if (surrogate == null) {
      res = new SurrogateGappedRegion(this, queryId, queryFrame, mSubjectFrames);
    } else {
      res =  ((SurrogateGappedRegion) surrogate).initialize(this, queryId, queryFrame);
    }
    mState = State.WRITTEN;
    return res;
  }

  @Override
  public void toString(final StringBuilder sb) {
    if (!isValid()) {
      sb.append("empty");
      return;
    }
    sb.append(mId);
    sb.append("[").append(mBuildStart).append("..").append(mBuildEnd).append("]");
    sb.append(mBuildSeqId);
    sb.append(":");
    sb.append("[").append(mQueryStart).append("..").append(mQueryEnd).append("]");
    sb.append("^").append(mBucket);
    sb.append("->").append(mNext.mId);
  }

  @Override
  public boolean integrity() {
    Exam.assertTrue(mWordSize > 0);
    if (mState == State.INIT || mState == State.WRITTEN) {
      Exam.assertTrue(mBuildSeqId >= 0);
      Exam.assertTrue(mBuildStart >= 0);
      //assertTrue(mBuildStart + mWordSize - 1 <= mBuildEnd);
      //assertTrue(mBuildStart % mStepSize == 0);
      //assertTrue((mBuildEnd - mWordSize + 1) % mStepSize == 0);
      Exam.assertTrue(mQueryStart >= 0);
      //assertTrue(mQueryStart <= mQueryEnd);
      Exam.assertTrue(mBuckets != null && mBucket == mBuckets.bucket(mBuildSeqId, mBuildEnd, mQueryEnd));
      Exam.assertTrue(mBucket >= 0);
      Exam.assertTrue(mNext != null);
    } else {
      Exam.assertEquals(NULL, mBuildSeqId);
      Exam.assertEquals(NULL, mBuildStart);
      Exam.assertEquals(NULL, mBuildEnd);
      Exam.assertEquals(NULL, mQueryStart);
      Exam.assertEquals(NULL, mQueryEnd);
      Exam.assertEquals(null, mBuckets);
      Exam.assertEquals(NULL, mBucket);
      Exam.assertEquals(null, mNext);
    }
    return true;
  }

}
