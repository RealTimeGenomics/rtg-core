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
package com.rtg.protein;

import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;

import com.rtg.reader.SequencesReader;
import com.rtg.util.diagnostic.Diagnostic;

/**
 * Groups length of reads into sizes that will be treated as equal for purposes of indexing
 */
public final class SequenceLengthBuckets {

  private static final int BUCKET_SIZE = 3;

  private final Map<Long, Long> mLengthCounts = new HashMap<>();
  private final long mMinLength;
  private final int mBucketSize;

  /**
   * Counts the number of sequences which fall into each protein length bucket
   * @param sr the sequence reader with mixed length sequences
   * @param minLength do not create buckets for lengths shorter than this.
   * @throws IOException when the sequence reader blows up.
   */
  public SequenceLengthBuckets(SequencesReader sr, long minLength) throws IOException {
    this(BUCKET_SIZE, minLength, sr.sequenceLengths(0, sr.numberSequences()));
  }

  private SequenceLengthBuckets(int bucketSize, long minLength, int[] lengths) {
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
