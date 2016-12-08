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

import com.rtg.launcher.BuildParams;
import com.rtg.reader.AlternatingSequencesReader;
import com.rtg.reader.SequencesReader;
import com.rtg.util.MathUtils;
import com.rtg.util.StringUtils;
import com.rtg.util.array.WrappedIntArray;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;

/**
 */
public final class GapBucketsInfo extends IntegralAbstract {

  private static final int INT_BYTES = 4;

  private final int mGap;

  private final int mStepSize;

  private final int mBucketScale;

  private final long[] mBucketLast;

  private final int mBits;

  private final int mLength;

  private final long mBucketMask;

  private final long mBucketTotal;

  private final boolean mReverseFrame;

  private final WrappedIntArray mReadLengths;

  /**
   * Build direct from <code>PositionParams</code>.
   * @param params all the necessary parameters.
   * @param scale how much to reduce the size of the buckets allocated.
   * @throws IOException If an I/O error occurs
   */
  public GapBucketsInfo(final PositionParams params, final int scale) throws IOException {
    this(params.build().sequences().reader(), params.build().stepSize(), params.build().sequences().mode().numberFrames(), scale, params.output().maxGap(), params.hashCountThreshold());
  }

  /**
   * Build direct from <code>PositionParams</code>.
   * @param params all the necessary parameters.
   * @param scale how much to reduce the size of the buckets allocated.
   * @param ngsMode always true, just for overloading
   * @throws IOException If an I/O error occurs
   */
  public GapBucketsInfo(final PositionParams params, final int scale, final boolean ngsMode) throws IOException {
    this(new AlternatingSequencesReader(params.ngsParams().buildFirstParams().reader(), params.ngsParams().buildSecondParams().reader()),
            params.build().stepSize(), params.build().sequences().mode().numberFrames(), scale, params.output().maxGap(), params.hashCountThreshold());
  }

  GapBucketsInfo(final BuildParams params, final int scale, final Integer gap, final Integer maxRepeat) throws IOException {
    this(params.sequences().reader(), params.stepSize(), params.sequences().mode().numberFrames(), scale, gap, maxRepeat);
  }

  //Only called directly from tests

  GapBucketsInfo(final SequencesReader reader, final int stepSize, final int numberFrames, final int scale, final Integer gap, final Integer maxRepeat) throws IOException {
    mReverseFrame = false;
    mReadLengths = null;
    mGap = gap;
    mBucketScale = scale;
    mStepSize = stepSize;
    //first work out how many buckets we need.
    final long start = 0;
    final long end = reader.numberSequences();
    Diagnostic.userLog("params for GapBuckets" + StringUtils.LS + "numSequences=" + reader.numberSequences() + " stepSize=" + stepSize + " numFrames=" + numberFrames);
    final long ns = reader.numberSequences()  * numberFrames;
    //System.err.println("ns=" + ns);
    if (ns > Integer.MAX_VALUE) {
      throw new IllegalArgumentException("too many sequences." + ns);
    }
    mBucketLast = new long[(int) ns];
    //a cumulative count of the sequence lengths will do it.
    long tot = 0;
    long max = 0;
    int j = 0;
    final int ss = mBucketScale * mStepSize;
    //System.err.println("bucketScale=" + mBucketScale + " stepSize=" + mStepSize + " ss=" + ss + " gap=" + mGap);
    //System.err.println("seqParams=" + seqParams);
    for (long i = start; i < end; ++i) {
      final int l = (reader.length(i) + mGap + ss - 1) / ss;
      //System.err.println("l=" + l + " length=" + reader.currentLength());
      for (int n = 0; n < numberFrames; ++n) {
        tot += l;
        if (l > max) {
          max = l;
        }
        mBucketLast[j] = tot - 1;
        ++j;
      }
    }

    mBucketTotal = tot;
    final long gappedLength = gap * (maxRepeat * 8L);
    final long targetLength = Math.min(Math.max(0, tot - 1), gappedLength);
    //make length of array a power of 2 so can do a mask rather than the much more expensive mod
    mBits = MathUtils.ceilPowerOf2Bits(targetLength);
    mLength = 1 << mBits;
    mBucketMask = mLength - 1;
    final String msg = "size of bucket array:" + mLength + " bits=" + mBits + " target length=" + targetLength + " gapped length=" + gappedLength + " tot=" + tot;
    //System.err.println(msg);
    Diagnostic.userLog(msg);
    Diagnostic.userLog("maximum subsequence:" + max);
    //System.err.println(this.toString());
    globalIntegrity();
  }

  /**
   * Get the actual number of buckets to be allocated.
   * This may be less than the bucket range because the range is deliberately
   * mapped to this smaller value.
   * @return the number of buckets to be allocated.
   */
  int length() {
    return mLength;
  }

  /**
   * Get the range of bucket values.
   * @return the range of bucket values.
   */
  long total() {
    return mBucketTotal;
  }

  /**
   * Get the mask used when accessing bucket indexes.
   * @return the 64 bit mask.
   */
  public long mask() {
    return mBucketMask;
  }

  /**
   * Get the first bucket assigned to the sequence.
   * @param seqid internal identifier of the sequence.
   * @return index of first assigned bucket.
   */
  long firstBucket(final int seqid) {
    if (seqid == 0) {
      return 0;
    }
    return mBucketLast[seqid - 1] + 1;
  }

  /**
   * Get the last bucket assigned to the sequence.
   * @param seqid internal identifier of the sequence.
   * @return index of last assigned bucket.
   */
  long lastBucket(final int seqid) {
    return mBucketLast[seqid];
  }

  /**
   * Get the bucket given the (template) sequence and a position in the
   * template and query sequence.
   * The bucket is fetched depending on the delta (difference between
   * template and query sequence positions).
   * @param seqid internal identifier of the sequence.
   * @param i position on the template sequence.
   * @param j  position on the query sequence.
   * @return index of the bucket.
   */
  long bucket(final int seqid, final int i, final int j) {
    final int delta;
    if (mReverseFrame) {
      delta = mReadLengths.get(seqid) - i - j;
    } else {
      delta = i - j;
    }
    //System.err.println("bucket delta=" + delta + " i=" + i + " j=" + j);
    final long b;
    final long first = firstBucket(seqid);
    final long last = lastBucket(seqid);
    assert first >= 0 && last >= first;
    final long length = last - first + 1;
    assert length >= 1;
    if (delta >= 0) {
      b = first + (delta / mStepSize) % length;
    } else {
      final int d = -(delta + 1) / mStepSize;
      assert d >= 0;
      b = last - (d % length);
    }
    assert b >= first && b <= last;
    //System.err.println("seqid=" + seqid + " i=" + i + " j=" + j + " delta=" + delta + " first=" + first + " last=" + last + " b=" + b);
    return b;
  }

  /**
   * Get the bucket given the (template) sequence and a position in the
   * template and query sequence.
   * The bucket is fetched depending on the delta (difference between
   * template and query sequence positions).
   * It differs from bucket in that it doesnt wrap at the ends so the
   * returned value may lie outside the first to last range.
   * @param seqid internal identifier of the sequence.
   * @param i position on the template sequence.
   * @param j  position on the query sequence.
   * @return index of the bucket.
   */
  long bucketEx(final int seqid, final int i, final int j) {
    final int delta = i - j;
    //System.err.println("bucket delta=" + delta + " i=" + i + " j=" + j);
    final long first = firstBucket(seqid);
    final long b;
    if (delta >= 0) {
      b = first + (delta / mStepSize);
    } else {
      final int d = -(delta + 1) / mStepSize;
      assert d >= 0;
      b = first - d - 1;
    }
    //System.err.println("seqid=" + seqid + " i=" + i + " j=" + j + " delta=" + delta + " first=" + first + " last=" + last + " b=" + b);
    return b;
  }

  /**
   * Get the total number of sequence ids that will be handled.
   * This may be more than the number of sequences in the reader because of framing.
   * Ony used for testing.
   * @return the total number of sequence ids that will be handled.
   */
  int numberSeqIds() {
    return mBucketLast.length;
  }

  /**
   *  Get number bytes used by this object.
   * @return number bytes used by this object.
   */
  public long bytes() {
    return mBucketLast.length * INT_BYTES;
  }

  /**
   * Whether buckets are running in reverse frame mode.
   * @return true if running in reverse
   */
  public boolean isReverseFrame() {
    return mReverseFrame;
  }

  @Override
  public void toString(final StringBuilder sb) {
    sb.append("GapBuckets pointers[").append(mBucketLast.length).append("]").append(StringUtils.LS);
    for (int i = 0; i < mBucketLast.length; ++i) {
      sb.append("[");
      sb.append(String.valueOf(i));
      sb.append("]->");
      sb.append(firstBucket(i));
      sb.append("..");
      sb.append(lastBucket(i));
      sb.append(StringUtils.LS);
      }
  }

  @Override
  public boolean globalIntegrity() {
    integrity();
    final int minBuckets = mGap / (mBucketScale * mStepSize);
    if (mBucketLast.length > 0) {
      long last = 0;
      for (int i = 1; i < mBucketLast.length; ++i) {
        final long l = mBucketLast[i];
        final long d = l - last;
        Exam.assertTrue(d >= minBuckets);
        last = l;
      }
    }
    return true;
  }

  @Override
  public boolean integrity() {
    Exam.assertTrue(0 < mStepSize);
    Exam.assertTrue(mGap > 0);
    Exam.assertTrue(mBucketLast != null);
    Exam.assertTrue(mLength > 0);
    Exam.assertTrue(mBits >= 0);
    Exam.assertTrue(mLength == 1L << mBits);
    final long mm = mLength - 1;
    Exam.assertTrue(((mLength & mBucketMask) == 0) && (mm & mBucketMask) == mm);
    return true;
  }
}
