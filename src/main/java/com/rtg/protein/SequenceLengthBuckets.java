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
package com.rtg.protein;

import java.util.HashMap;
import java.util.Map;
import java.util.Set;

import com.rtg.util.diagnostic.Diagnostic;

/**
 * Groups length of reads into sizes that will be treated as equal for purposes of indexing
 */
public final class SequenceLengthBuckets {

  private final Map<Long, Long> mLengthCounts = new HashMap<>();
  private final long mMinLength;
  private final int mBucketSize;

  /**
   * Counts the number of sequences which fall into each protein length bucket
   * @param bucketSize the size of each bucket
   * @param minLength do not create buckets for lengths shorter than this
   * @param lengths the lengths to be bucketed
   */
  SequenceLengthBuckets(int bucketSize, long minLength, int[] lengths) {
    mBucketSize = bucketSize;
    mMinLength = minLength;
    for (long length : lengths) {
      if (length >= minLength) {
        final long bucket = bucket(length);
        // System.err.println("length: " + length + " bucket: " + bucket);
        if (mLengthCounts.containsKey(bucket)) {
          mLengthCounts.put(bucket, mLengthCounts.get(bucket) + 1);
        } else {
          mLengthCounts.put(bucket, 1L);
        }
      }
    }
    Diagnostic.developerLog(this.toString());
  }

  /**
   * Determine the bucket a read of a specific length will fall into
   * @param length the length of a read
   * @return the bucket id or -1 if the read is shorter than the minimum length specified at construction time.
   */
  public long bucket(long length) {
    if (length < mMinLength) {
      return -1;
    }
    return (length / mBucketSize) * mBucketSize;
  }

  /**
   *
   * @return the set of read lengths for which buckets exist
   */
  public Set<Long> getBuckets() {
    return mLengthCounts.keySet();
  }

  /**
   * @param bucketId which bucket would you like to count?
   * @return the number of reads which fall into the bucket
   */
  public long getCount(long bucketId) {
    if (!mLengthCounts.containsKey(bucketId)) {
      return 0;
    }
    return mLengthCounts.get(bucketId);
  }

  @Override
  public String toString() {
    return "SequenceLengthBuckets: bucketSize=" + mBucketSize + " minLength=" + mMinLength + " buckets=" + mLengthCounts;
  }
}
