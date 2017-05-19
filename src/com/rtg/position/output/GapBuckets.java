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

import com.rtg.util.StringUtils;
import com.rtg.util.array.objectindex.ObjectCreate;
import com.rtg.util.array.objectindex.ObjectIndex;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;


/**
 * Holds bucket information for GappedOutput.
 * @param <A> type of stored objects.
 */
public class GapBuckets<A> extends IntegralAbstract {

  private static final int OBJECT_BYTES = 4;

  private final GapBucketsInfo mInfo;

  private final long mBucketMask;

  private final ObjectIndex<A> mBuckets;

  /**
   * @param info information needed for construction.
   */
  GapBuckets(final GapBucketsInfo info) {
    mInfo = info;
    mBuckets = ObjectCreate.createIndex(info.length());
    mBucketMask = info.mask();
    //System.err.println(this.toString());
    assert globalIntegrity();
    assert isEmpty();
  }

  /**
   * Get the pointer for the specified bucket.
   * @param bucket index.
   * @return the pointer.
   */
  A get(final long bucket) {
    return mBuckets.get(bucket & mBucketMask);
  }

  /**
   * Set the value of the specified bucket.
   * @param bucket index.
   * @param region to be stored in bucket (may be null).
   */
  void set(final long bucket, final A region) {
    mBuckets.set(bucket & mBucketMask, region);
  }

  /**
   * Get the first bucket assigned to the sequence.
   * @param seqId internal identifier of the sequence.
   * @return index of first assigned bucket.
   */
  long firstBucket(final int seqId) {
    return mInfo.firstBucket(seqId);
  }

  /**
   * Get the last bucket assigned to the sequence.
   * @param seqId internal identifier of the sequence.
   * @return index of last assigned bucket.
   */
  long lastBucket(final int seqId) {
    return mInfo.lastBucket(seqId);
  }

  /**
   * Get the bucket given the (template) sequence and a position in the
   * template and query sequence.
   * The bucket is fetched depending on the delta (difference between
   * template and query sequence positions).
   * @param seqId internal identifier of the sequence.
   * @param i position on the template sequence.
   * @param j  position on the query sequence.
   * @return index of the bucket.
   */
  long bucket(final int seqId, final int i, final int j) {
    return mInfo.bucket(seqId, i, j);
  }

  /**
   * Get the bucket given the (template) sequence and a position in the
   * template and query sequence.
   * The bucket is fetched depending on the delta (difference between
   * template and query sequence positions).
   * It differs from bucket in that it doesnt wrap at the ends so the
   * returned value may lie outside the first to last range.
   * @param seqId internal identifier of the sequence.
   * @param i position on the template sequence.
   * @param j  position on the query sequence.
   * @return index of the bucket.
   */
  long bucketEx(final int seqId, final int i, final int j) {
    return mInfo.bucketEx(seqId, i, j);
  }

  /**
   * Check if all the buckets are empty.
   * @return true iff all the buckets are empty.
   */
  boolean isEmpty() {
    for (long i = 0; i < mBuckets.length(); ++i) {
      final A region = mBuckets.get(i);
      if (region != null) {
        return false;
      }
    }
    return true;
  }

  /**
   * Get the total number of sequence ids that will be handled.
   * This may be more than the number of sequences in the reader because of framing.
   * Ony used for testing.
   * @return the total number of sequence ids that will be handled.
   */
  int numberSeqIds() {
    return mInfo.numberSeqIds();
  }

  /**
   * Get the total number of buckets.
   * Only used for integrity checks.
   * @return the total number of buckets.
   */
  long numberBuckets() {
    return mBuckets.length();
  }

  /**
   *  Get number bytes used by this object.
   * @return number bytes used by this object.
   */
  public long bytes() {
    return mBuckets.length() * OBJECT_BYTES;
  }

  @Override
  public void toString(final StringBuilder sb) {
    mInfo.toString(sb);
    sb.append("buckets [").append(mBuckets.length()).append("]").append(StringUtils.LS);
    for (long i = 0; i < mBuckets.length(); ++i) {
      final A a = mBuckets.get(i);
      if (a != null) {
        sb.append("[").append(i).append("] ").append(a.toString()).append(StringUtils.LS);
      }
    }
  }

  @Override
  public final boolean globalIntegrity() {
    integrity();
    mInfo.globalIntegrity();
    return true;
  }

  @Override
  public boolean integrity() {
    Exam.assertNotNull(mBuckets);
    if (mBuckets != null) {
      if (mInfo.numberSeqIds() == 0) {
        Exam.assertTrue("" + mBuckets.length(), mBuckets.length() == 1);
      } else {
        final boolean check = mInfo.lastBucket(numberSeqIds() - 1) == mInfo.total() - 1;
        Exam.assertTrue(this.toString(), check);
      }
    }
    return true;
  }
}
