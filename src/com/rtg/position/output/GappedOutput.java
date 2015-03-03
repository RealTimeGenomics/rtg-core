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

package com.rtg.position.output;

import java.io.IOException;
import java.util.HashSet;
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

  private final AbstractGappedRegion<G>[][] mRegions;

  private final int mNumberColumns;

  private final int[] mNextRegion;

  private final GappedRegionFactory<G> mRegionFactory;

  private final PositionWriter mWriter;

  private int mCurrColumn = NULL;

  private int mId = 0;

  private final boolean mReverseFrame;

  private final int mRepeatThreshold;

  //hang on to these for clone calls

  private final GapBucketsInfo mBucketInfo;

  /**
   *
   * @param params parameters.
   * @param regionFactory factory for creating regions
   * @param scorer scoring implementation used
   * @param buildLengths lengths of sequences.
   * @param out where to write the output.
   * @param writer where to write output for unmapped reads
   * @param bucketInfo bucket info
   * @param repeatFreq the repeat frequency threshold
   */
  GappedOutput(final PositionParams params, GappedRegionFactory<G> regionFactory, GapScorer scorer, final ImmutableIntArray buildLengths, final Appendable out, final PositionWriter writer, final GapBucketsInfo bucketInfo, int repeatFreq) {
    super(params, out);
    mRepeatThreshold = repeatFreq;
    mBucketInfo = bucketInfo;
    mBuildLengths = buildLengths;
    final Integer maxGap = params.output().maxGap();
    //System.err.println("params=" + params);
    //System.err.println(maxGap + "=maxGap");
    mMaxGap = maxGap;
    final BuildParams buildParams = params.build();
    mWordSize = buildParams.windowSize();
    mProbabilities = scorer; // params.output().format().getScorer(params);
    //System.err.println("probabilities=" + mProbabilities + " GappedOutput");
    mMinDelta = mProbabilities.minDelta();
    mMaxDelta = mProbabilities.maxDelta();
    Diagnostic.developerLog("GappedOutput minDelta=" + mMinDelta + " maxDelta=" + mMaxDelta + " maxGap=" + maxGap);
    mReverseFrame = bucketInfo.isReverseFrame();
    mBuckets = new GapBuckets<>(bucketInfo);
    mNumberColumns = maxGap + 1;
    mRegions = regionArray();
    mRegionFactory = regionFactory;
    mNextRegion = new int[mNumberColumns];
    mWriter = writer;
    mScanner = new GappedOutputScannerImpl<>(this);
  }

  private AbstractGappedRegion<G>[][] regionArray() {
    //System.err.println(" threshold=" + params.threshold());
    @SuppressWarnings("unchecked")
    final AbstractGappedRegion<G>[][] ret = (AbstractGappedRegion<G>[][]) new AbstractGappedRegion<?>[mNumberColumns][];
    for (int i = 0; i < ret.length; i++) {
      @SuppressWarnings("unchecked")
      final AbstractGappedRegion<G>[] gs = (AbstractGappedRegion<G>[]) new AbstractGappedRegion<?>[mRepeatThreshold];
      ret[i] = gs;
    }
    return ret;
  }

  @Override
  public void endQuery() throws IOException {
    //System.err.println("endQuery");
    //printRegions(System.err);
    if (mSearchPosition != NULL) {
      for (int i = mSearchPosition + 1; i <= mSearchPosition + mNumberColumns; i++) {
        final int im = i % mNumberColumns;
        //System.err.println("flushing i=" + i + " im=" + im);
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
    for (int i = 0; i < mRegions.length; i++) {
      for (int j = 0; j < mRegions[i].length; j++) {
        assert mRegions[i][j] == null || !mRegions[i][j].isValid() : "endQuery i=" + i + " j=" + j + " region=" + mRegions[i][j];
      }
    }
    return true;
  }

//  void printRegions(final PrintStream out) {
//    out.println("Regions[" + mRegions.length + "][" + mRegions[0].length + "] search=" + mSearchPosition + " curr. col=" + mCurrColumn);
//    for (int i = 0; i < mRegions.length; i++) {
//      //System.err.println("i=" + i);
//      out.print("[" + i + "] ");
//      for (int j = 0; j < mRegions[i].length; j++) {
//        if (mRegions[i][j].isValid()) {
//          out.print((j == 0 ? "" : ",") + " [" + j + "]" + mRegions[i][j]);
//        }
//      }
//      out.println();
//    }
//  }

  @Override
  public void hit(final int seqId, final int posn) throws IOException {
    final int index = mNextRegion[mCurrColumn];
    //System.err.println("hit seqId=" + seqId + " posn=" + posn + " curr=" + mCurrColumn + " index=" + index);
    @SuppressWarnings("unchecked")
    final G cu = (G) mRegions[mCurrColumn][index];
    final G curr;
    if (cu == null) { //lazily contruct regions
      mId++;
      //System.err.println("buildLengths=" + mBuildLengths);
      curr = mRegionFactory.region(mId, mProbabilities, mParams, mBuildLengths);
      mRegions[mCurrColumn][index] = curr;
    } else {
      curr = cu;
    }
    final long bucket = mBuckets.bucket(seqId, posn, mSearchPosition);
    //System.err.println("bucket=" + bucket + " seqId=" + seqId);
    final G last = mBuckets.get(bucket);
    final G next;
    if (last == null) {
      next = curr;
    } else {
      next = last.next();
      last.setNext(curr);
    }
    curr.initialize(seqId, posn, mSearchPosition, mBuckets, next, mReverseFrame);
    mNextRegion[mCurrColumn]++;
    mBuckets.set(bucket, curr);
    super.hit(seqId, posn);
  }

  @Override
  public void setPosition(final int position) throws IOException {
    //System.err.println("setPosition position=" + position);
    if (mSearchPosition != NULL) {
      if (position == mSearchPosition + 1) {
        final int in = mCurrColumn + 1;
        final int im = in == mNumberColumns ? 0 : in;
        flush(im);
        mCurrColumn = im;
      } else {
        final int end = mSearchPosition + mNumberColumns;
        final int hi = position >= end ? end : position;
        //System.err.println("mSearchPosition=" + mSearchPosition + " hi=" + hi);
        for (int i = mSearchPosition + 1; i <= hi; i++) {
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
    //if (mNextRegion[col] > 0) {
    //System.err.println("flush col=" + col + " next=" + mNextRegion[col]);
    //}
    //System.err.println(this.toString());
    for (int i = 0; i < mNextRegion[col]; i++) {
      //System.err.println("flush col=" + col + " i=" + i);
      final AbstractGappedRegion<G> region = mRegions[col][i];
      //this is tricky
      //the bucket may not point to the region we have just found
      //so the merged region may be another one - MUST keep going here
      //so we dont miss because the other is no longer initialized by the time we get to it
      final long b = region.bucket();
      //System.err.println(b + "=b");
      while (region.isInitialized()) {
        final G actual = mBuckets.get(b).next();
        flushBucket(actual);
        assert actual.isValid() && !actual.isInitialized();
        removeBucket(b);
        assert actual.isValid() && !actual.isInitialized();
        actual.reset();
      }
    }
    mNextRegion[col] = 0;
    //System.err.println(this.toString());
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
      //System.err.println("region=" + region + " best=" + best + " bestBest=" + bestBest);
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
    return mBuckets.bytes() + mRegions.length * mRegions[0].length * 42L + mNextRegion.length * 4;
  }

  @Override
  public void memToString(final StringBuilder sb) {


  }

  @Override
  public void toString(final StringBuilder sb) {
    sb.append("GappedOutput wordSize=").append(mWordSize).append(" maxGap=").append(mMaxGap).append(" minDelta=").append(mMinDelta).append(" maxDelta=").append(mMaxDelta).append(" score=").append(score()).append(StringUtils.LS);
    for (int i = 0; i < mRegions.length; i++) {
      sb.append("[").append(i).append("]").append(mNextRegion[i]);
      for (int j = 0; j < mNextRegion[i]; j++) {
        sb.append(" ").append(mRegions[i][j]);
      }
      sb.append(StringUtils.LS);
    }
  }

  @Override
  public boolean globalIntegrity() {
    integrity();
    //scan the regions array
    if (!mReverseFrame) {
      for (int i = 0; i < mRegions.length; i++) {
        final AbstractGappedRegion<G>[] column = mRegions[i];
        final int init = mNextRegion[i];
        Exam.assertTrue(init >= 0 && init <= column.length);
        int qEnd = -1;
        for (int j = 0; j < init; j++) {
          final AbstractGappedRegion<G> region = column[j];
          region.integrity();
          Exam.assertTrue(region.isInitialized());
          final int k = (region.queryEnd() - mWordSize + 1) % (mMaxGap + 1);
          //System.err.println("i=" + i + " k=" + k + " region=" + region);
          Exam.assertEquals(i, k);
          if (qEnd == -1) {
            qEnd = region.queryEnd();
          } else {
            Exam.assertEquals(qEnd, region.queryEnd());
          }
          final long bucket = region.bucket();
          Exam.assertTrue(mBuckets.get(bucket) != null);
        }
        for (int j = init; j < column.length; j++) {
          final AbstractGappedRegion<G> region = column[j];
          if (region != null) {
            region.integrity();
            Exam.assertTrue(!region.isInitialized());
          }
        }
      }

      //scan all the buckets
      for (int i = 0; i < mBuckets.numberBuckets(); i++) {
        final AbstractGappedRegion<G> last = mBuckets.get(i);
        if (last == null) {
          continue;
        }
        final Set<Integer> ids = new HashSet<>();
        G ptr = last.next();
        int lastEnd = 0;
        while (true) {
          //all point to same bucket
          //System.err.println("i=" + i + " ptr.bucket()=" + ptr.bucket() + " mod=" + (ptr.bucket() % mBuckets.numberBuckets()));
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

    public GappedOutputScannerImpl(final GappedOutput<G> outer) {
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
      //System.err.println("scan region=" + region + " forward=" + forward + " lo=" + lo + " hi=" + hi);
      for (long b = lo; b <= hi; b++) {
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
            //System.err.println(" score=" + score + " ptr=" + ptr);
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
    return new GappedOutput<>(mParams, mRegionFactory, mProbabilities, mBuildLengths, mOut, mWriter.reverseClone(), mBucketInfo, mRepeatThreshold);
  }

}
