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

import com.rtg.util.StringUtils;
import com.rtg.util.array.objectindex.ObjectCreate;
import com.rtg.util.array.objectindex.ObjectIndex;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;


/**
 * Holds bucket information for GappedOutput.
 * @param <A> type of stored objects.
 */
public final class GapBuckets<A> extends IntegralAbstract {

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
        sb.append("[").append(i).append("] ").append(a).append(StringUtils.LS);
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
        Exam.assertTrue("length=" + mBuckets.length(), mBuckets.length() == 1);
      } else {
        final boolean check = mInfo.lastBucket(numberSeqIds() - 1) == mInfo.total() - 1;
        Exam.assertTrue(this.toString(), check);
      }
    }
    return true;
  }
}
