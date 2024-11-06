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

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import com.rtg.launcher.BuildParams;
import com.rtg.mode.Frame;
import com.rtg.util.StringUtils;
import com.rtg.util.array.ImmutableIntArray;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.integrity.Exam;

/**
 * @param <G> type of gapped HashingRegion
 */
public class GappedOutput <G extends AbstractGappedRegion<G>> extends AbstractPositionOutput {

  private final ImmutableIntArray mBuildLengths;

  private final int mWordSize;

  private final Integer mMaxGap;

  private final GapScorer mProbabilities;

  private final int mMinDelta;

  private final int mMaxDelta;

  private final GapBuckets<G> mBuckets;

  private final List<AbstractGappedRegion<G>>[] mRegions;

  private final int mNumberColumns;

  private final GappedRegionFactory<G> mRegionFactory;

  private final PositionWriter mWriter;

  private int mCurrColumn = NULL;

  private int mId = 0;

  private final boolean mReverseFrame;

  //hang on to these for clone calls

  private final GapBucketsInfo mBucketInfo;

  /**
   * @param params parameters.
   * @param regionFactory factory for creating regions
   * @param scorer scoring implementation used
   * @param buildLengths lengths of sequences.
   * @param out where to write the output.
   * @param writer where to write output for unmapped reads
   * @param bucketInfo bucket info
   */
  GappedOutput(final PositionParams params, GappedRegionFactory<G> regionFactory, GapScorer scorer, final ImmutableIntArray buildLengths, final Appendable out, final PositionWriter writer, final GapBucketsInfo bucketInfo) {
    super(params, out);
    mBucketInfo = bucketInfo;
    mBuildLengths = buildLengths;
    final Integer maxGap = params.output().maxGap();
    mMaxGap = maxGap;
    final BuildParams buildParams = params.build();
    mWordSize = buildParams.windowSize();
    mProbabilities = scorer;
    mMinDelta = mProbabilities.minDelta();
    mMaxDelta = mProbabilities.maxDelta();
    Diagnostic.developerLog("GappedOutput minDelta=" + mMinDelta + " maxDelta=" + mMaxDelta + " maxGap=" + maxGap);
    mReverseFrame = bucketInfo.isReverseFrame();
    mBuckets = new GapBuckets<>(bucketInfo);
    mNumberColumns = maxGap + 1;
    mRegions = regionArray();
    mRegionFactory = regionFactory;
    mWriter = writer;
    mScanner = new GappedOutputScannerImpl<>(this);
  }

  private List<AbstractGappedRegion<G>>[] regionArray() {
    @SuppressWarnings("unchecked")
    final List<AbstractGappedRegion<G>>[] ret = (List<AbstractGappedRegion<G>>[]) new List<?>[mNumberColumns];
    for (int i = 0; i < ret.length; ++i) {
      ret[i] = new ArrayList<>();
    }
    return ret;
  }

  @Override
  public void endQuery() throws IOException {
    if (mSearchPosition != NULL) {
      for (int i = mSearchPosition + 1; i <= mSearchPosition + mNumberColumns; ++i) {
        final int im = i % mNumberColumns;
        flush(im);
      }
    }
    mCurrColumn = NULL;
    assert checkRegionsEmpty();
    super.endQuery();
  }

  @Override
  public void endQuerySequence() throws IOException {
    mWriter.endQuery(mSearchSeqId);
    super.endQuerySequence();
  }


  @Override
  public void endAll() throws IOException {
    super.endAll();
    mWriter.endAll();
    //TODO write score
    mProbabilities.close();
  }

  private boolean checkRegionsEmpty() {
    for (int i = 0; i < mRegions.length; ++i) {
      for (int j = 0; j < mRegions[i].size(); ++j) {
        assert mRegions[i].get(j).isValid() : "endQuery i=" + i + " j=" + j + " region=" + mRegions[i].get(j);
      }
    }
    return true;
  }

  @Override
  public void hit(final int seqId, final int posn) throws IOException {
    ++mId;
    final G curr = mRegionFactory.region(mId, mProbabilities, mParams, mBuildLengths);
    mRegions[mCurrColumn].add(curr);
    final long bucket = mBuckets.bucket(seqId, posn, mSearchPosition);
    final G last = mBuckets.get(bucket);
    final G next;
    if (last == null) {
      next = curr;
    } else {
      next = last.next();
      last.setNext(curr);
    }
    curr.initialize(seqId, posn, mSearchPosition, mBuckets, next, mReverseFrame);
    mBuckets.set(bucket, curr);
    super.hit(seqId, posn);
  }

  @Override
  public void setPosition(final int position) throws IOException {
    if (mSearchPosition != NULL) {
      if (position == mSearchPosition + 1) {
        final int in = mCurrColumn + 1;
        final int im = in == mNumberColumns ? 0 : in;
        flush(im);
        mCurrColumn = im;
      } else {
        final int end = mSearchPosition + mNumberColumns;
        final int hi = position >= end ? end : position;
        for (int i = mSearchPosition + 1; i <= hi; ++i) {
          final int im = i % mNumberColumns;
          flush(im);
        }
        mCurrColumn = position % mNumberColumns;
      }
    } else {
      mCurrColumn = position % mNumberColumns;
    }
    assert mState == State.POSITION;
    mSearchPosition = position;
    mState = State.HIT;
  }

  private void flush(final int col) throws IOException {
    for (int i = 0; i < mRegions[col].size(); ++i) {
      final AbstractGappedRegion<G> region = mRegions[col].get(i);
      //this is tricky
      //the bucket may not point to the region we have just found
      //so the merged region may be another one - MUST keep going here
      //so we dont miss because the other is no longer initialized by the time we get to it
      final long b = region.bucket();
      while (region.isInitialized()) {
        final G actual = mBuckets.get(b).next();
        flushBucket(actual);
        assert actual.isValid() && !actual.isInitialized();
        removeBucket(b);
        assert actual.isValid() && !actual.isInitialized();
        actual.reset();
      }
    }
    mRegions[col].clear();
    assert globalIntegrity();
  }

  private void flushBucket(final G region) throws IOException {
    scanAll(region, true);
    final G best = mScanner.result();
    final double bestScore = mScanner.resultScore();
    if (bestScore != 0.0) {
      final G bestBest;
      if (best != null) {
        scanAll(best, false);
        bestBest =  mScanner.result();
      } else {
        bestBest = null;
      }
      if (best == null || region != bestBest) { //no reciprocal match write out the region
        mWriter.write(region, mSearchSeqId, mSearchFrame, mQueryLength, mQueryEffectiveLength);
        assert region.isValid() && !region.isInitialized();
        return;
      }
    }
    best.merge(region);
    region.write(); //set state
  }

  /*
   * Remove the first region from the front of the bucket.
   * Carefully deals with the case when this is the last region in the bucket chain.
   */
  private void removeBucket(final long bucket) {
    final G last = mBuckets.get(bucket);
    assert last != null;
    final G first = last.next();
    if (last == first) { //remove the last one
      mBuckets.set(bucket, null);
    } else {
      final G second = first.next();
      last.setNext(second);
    }
  }

  void scanAll(final G region, final boolean forward) {
    final long lo = region.bucketEx(mMinDelta);
    final long hi = region.bucketEx(mMaxDelta);
    //System.err.println("bestMatch mMinDelta=" + mMinDelta + " mMaxDelta=" + mMaxDelta + " region=" + region);
    final long first = region.firstBucket();
    final long last = region.lastBucket();
    mScanner.scanAll(region, forward, first, last, lo, hi);
    //ignore return value since we need more than one from it anyway
  }

  @Override
  public double score() {
    return mWriter.score();
  }

  final Scanner<G> mScanner;

  @Override
  public long bytes() {
    long total = mBuckets.bytes();
    for (final List<AbstractGappedRegion<G>> mRegion : mRegions) {
      total += mRegion.size() * 42L + 4;
    }
    return total;
  }

  @Override
  public void memToString(final StringBuilder sb) {


  }

  @Override
  public void toString(final StringBuilder sb) {
    sb.append("GappedOutput wordSize=").append(mWordSize).append(" maxGap=").append(mMaxGap).append(" minDelta=").append(mMinDelta).append(" maxDelta=").append(mMaxDelta).append(" score=").append(score()).append(StringUtils.LS);
    for (int i = 0; i < mRegions.length; ++i) {
      sb.append("[").append(i).append("]").append(mRegions[i].size());
      for (int j = 0; j < mRegions[i].size(); ++j) {
        sb.append(" ").append(mRegions[i].get(j));
      }
      sb.append(StringUtils.LS);
    }
  }

  @Override
  public boolean globalIntegrity() {
    integrity();
    //scan the regions array
    if (!mReverseFrame) {
      for (int i = 0; i < mRegions.length; ++i) {
        final List<AbstractGappedRegion<G>> column = mRegions[i];
        final int init = column.size();
        int qEnd = -1;
        for (final AbstractGappedRegion<G> region : column) {
          region.integrity();
          Exam.assertTrue(region.isInitialized());
          final int k = (region.queryEnd() - mWordSize + 1) % (mMaxGap + 1);
          Exam.assertEquals(i, k);
          if (qEnd == -1) {
            qEnd = region.queryEnd();
          } else {
            Exam.assertEquals(qEnd, region.queryEnd());
          }
          final long bucket = region.bucket();
          Exam.assertTrue(mBuckets.get(bucket) != null);
        }
        for (int j = init; j < column.size(); ++j) {
          final AbstractGappedRegion<G> region = column.get(j);
          assert region != null;
          region.integrity();
          Exam.assertTrue(!region.isInitialized());
        }
      }

      //scan all the buckets
      final Set<Integer> ids = new HashSet<>();
      for (int i = 0; i < mBuckets.numberBuckets(); ++i) {
        final AbstractGappedRegion<G> last = mBuckets.get(i);
        if (last == null) {
          continue;
        }
        ids.clear();
        G ptr = last.next();
        int lastEnd = 0;
        while (true) {
          //all point to same bucket
          Exam.assertEquals(i, ptr.bucket() % mBuckets.numberBuckets());
          //end point must increase along chain
          Exam.assertTrue(ptr.queryEnd() >= lastEnd);
          lastEnd = ptr.queryEnd();
          Exam.assertTrue(ids.add(ptr.id()));
          if (ptr == last) {
            break;
          }
          ptr = ptr.next();
        }
      }
    }
    return true;
  }

  @Override
  public boolean integrity() {
    super.integrity();
    if (mState == State.SEQUENCE) {
      //nothing
    } else if (mState == State.QUERY) {
      Exam.assertEquals(NULL, mCurrColumn);
    } else {
      Exam.assertTrue(mCurrColumn >= 0 && mCurrColumn < mRegions.length);
    }
    if (mSearchPosition != NULL) {
      Exam.assertEquals(mCurrColumn, mSearchPosition % mNumberColumns);
    }
    return true;
  }



  private static class GappedOutputScannerImpl<G extends AbstractGappedRegion<G>> extends Scanner<G> {

    GappedOutput<G> mOuter;
    private double mBestScore = Double.NEGATIVE_INFINITY;
    G mBestRegion = null;

    GappedOutputScannerImpl(final GappedOutput<G> outer) {
      super();
      this.mOuter = outer;
    }

    @Override
    protected void init() {
      mBestScore = Double.NEGATIVE_INFINITY;
      mBestRegion = null;
    }

    @Override
    protected G result() {
      return mBestRegion;
    }

    @Override
    protected double resultScore() {
      return mBestScore;
    }



    @Override
    protected void scan(final G region, final boolean forward, final long lo, final long hi) {
      for (long b = lo; b <= hi; ++b) {
        //have a bucket - scan all the regions linked to the bucket
        final G lastRegion = mOuter.mBuckets.get(b);
        if (lastRegion != null) {
          final G firstRegion = lastRegion.next();
          assert firstRegion != null;
          G ptr = firstRegion;
          while (true) {
            final double score;
            if (forward) {
              score = ptr.gapScore(mOuter.mProbabilities, region);
            } else {
              if (ptr.queryEnd() <= region.queryEnd()) {
                score = region.gapScore(mOuter.mProbabilities, ptr);
              } else {
                score = Double.NEGATIVE_INFINITY;
              }
            }
            if (score > mBestScore) {
              mBestScore = score;
              mBestRegion = ptr;
            }
            if (mBestScore == 0.0 || ptr == lastRegion) {
              break;
            }
            ptr = ptr.next();
          }
        }
      }
    }
  }

  @Override
  public void nextQuery(final Frame frame, final int seqId) {
    super.nextQuery(frame, seqId);
    mProbabilities.setQueryId(seqId, frame);
  }

  @Override
  public PositionOutput reverseClone() {
    return new GappedOutput<>(mParams, mRegionFactory, mProbabilities, mBuildLengths, mOut, mWriter.reverseClone(), mBucketInfo);
  }

}
