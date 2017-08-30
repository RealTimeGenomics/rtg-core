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
package com.rtg.index;

import java.io.IOException;
import java.util.Arrays;

import com.rtg.index.params.CreateParams;
import com.rtg.util.array.ArrayUtils;
import com.rtg.util.array.IndexSorter;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.OneShotTimer;
import com.rtg.util.integrity.Exam;

/**
 * Does not compress the hash codes.
 *
 */
public class IndexSimple extends IndexBase {

  //Various statistic counters set during searching

  /** Total number of attempts to find a hash code. */
  private long mSearchCount = 0;

  /** Number of searches that see a miss in the bit vector table. */
  private long mBitVectorMissCount = 0;

  /** Number of searches that see a hit in the bit vector table. */
  private long mBitVectorHitCount = 0;

  /** Number of searches that get a miss in the initial position table. */
  private long mMissCount = 0;

  /** Number of searches that get a hit in the initial position table. */
  private long mInitialHitCount = 0;

  /** Number of searches that get past frequency threshold. */
  private long mFinalHitCount = 0;

  /** Number of results returned as hits. */
  private long mTotalFinds = 0;

  /**
   * Constructs an empty index.
   *
   * @param indexParams holds all the values needed for constructing the index.
   * @param filter the hash filter
   * @param numberThreads number of threads appropriate for parallel execution.
   */
  public IndexSimple(final CreateParams indexParams, IndexFilterMethod filter, final int numberThreads) {
    super(indexParams, filter, numberThreads);
    if (mHashBits < 1 || mHashBits > 64) {
      throw new RuntimeException("Hash bits set to invalid value bits=" + mHashBits);
    }

    assert !mParams.compressHashes();
  }

  @Override
  public final void add(final long hash, final long value) {
    //System.err.println(mState + " add(" + hash + ", " + value + ")");

    assert mState == IndexState.ADD;
    assert value >= 0;
    if (mInitialHashes >= mHash.length()) {
      throw new RuntimeException("Too many items added:" + mInitialHashes);
    } else {
      mHash.set(mInitialHashes, hash);
      mValue.set(mInitialHashes, value);
    }
    ++mInitialHashes;
  }

  @Override
  protected final long position(final long hash) {
    final long k = (hash >> mSRBits) + mIncrement;
    assert k >= 0 && k < mInitialPositionLength - 1 : "k=" + k + " hash=" + hash + " srbits=" + mSRBits;
    return k;
  }

  @Override
  public void freeze() {
    assert mState == IndexState.ADD;
    final OneShotTimer sort = new OneShotTimer("Index_sort");
    IndexSorter.sort(mHash, mValue, mInitialHashes);
    sort.stopLog();
    //fill overflow table
    final OneShotTimer over = new OneShotTimer("Index_pointer");
    createPointers();
    over.stopLog();

    final OneShotTimer posn = new OneShotTimer("Index_position");
    populateInitialPointer();
    posn.stopLog();
    if (mHashVector != null) {
      final OneShotTimer bitv = new OneShotTimer("Index_bitVector");
      //set bit vector
      for (long i = 0; i < mNumValues; ++i) {
        final long hash = mHash.get(i);
        mHashVector.set(hash);
      }
      bitv.stopLog();
    }

    checkRepeatFrequency();

    mState = IndexState.FROZEN;
    assert globalIntegrity();
  }

  private void populateInitialPointer() {
    long lastKey = 0;
    long count = 0;
    if (mNumHashes == 0) {
      return;
    }

    for (long i = 0; i < mNumValues;) {
      final long hash = mHash.get(i);
      final long key = position(hash);
      assert key >= 0 && key < mInitialPositionLength - 1;
      if (key > lastKey && lastKey > 0) {
        //System.err.println("bucket 1 count=" + count);
        ++mBucketCount0;
        assert count >= 1; // : "count=" + count + " i=" + i;
        mBucketCount1 += count;
        mBucketCount2 += count * count;
        count = 0;
      }
      for (long j = lastKey + 1; j <= key; ++j) {
        mInitialPosition.set(j, i);
      }
      lastKey = key;
      long l = 0;
      while (i < mNumValues && hash == mHash.get(i)) {
        ++l;
        ++i;
      }
      count += l;
    }
    if (count > 0) {
      //System.err.println("bucket 2 count=" + count);
      ++mBucketCount0;
      assert count >= 1 : count;
      mBucketCount1 += count;
      mBucketCount2 += count * count;
    }
    for (long j = lastKey + 1; j < mInitialPositionLength; ++j) {
      mInitialPosition.set(j, mNumValues);
    }
  }

  @Override
  public SparseFrequencyHistogram getSparseFrequencyHistogram() {
    int[] freqDist = new int[1024];
    SparseFrequencyHistogram freqHist = new SparseFrequencyHistogram();
    long mergeCount = 0;
    int numUsed = 0;
    for (long i = 0; i < mInitialHashes;) {
      final long hash = mHash.get(i);
      int freq = 1;
      ++i;
      while (i < mInitialHashes && hash == mHash.get(i)) {
        ++i;
        ++freq;
      }
      if (numUsed >= freqDist.length) {
        if (freqDist.length == MAX_FREQ_DIST_SIZE) {
          freqHist = SparseFrequencyHistogram.merge(freqHist, SparseFrequencyHistogram.fromIndividualFrequencies(freqDist, numUsed));
          ++mergeCount;
          numUsed = 0;
        } else {
          int length = freqDist.length * 3 / 2;
          if (length < 0 || length > MAX_FREQ_DIST_SIZE) {
            length = MAX_FREQ_DIST_SIZE;
          }
          freqDist = Arrays.copyOf(freqDist, length);
        }
      }
      freqDist[numUsed++] = freq;
    }
    freqHist = SparseFrequencyHistogram.merge(freqHist, SparseFrequencyHistogram.fromIndividualFrequencies(freqDist, numUsed));
    ++mergeCount;
    Diagnostic.developerLog("IndexSimple " + mergeCount + " Frequency Histogram Merges");
    return freqHist;
  }

  @Override
  protected void compact() {
    mMaxHashCount = 0;

    mIndexFilterMethod.initialize(this);
    long prevk = 0;
    for (long i = 0; i < mInitialHashes;) {
      long k = prevk;
      final long hash = mHash.get(i);
      while (i < mInitialHashes && hash == mHash.get(i)) {
        mHash.set(k, hash);
        mValue.set(k, mValue.get(i));
        ++k;
        ++i;
      }
      final long nh = k - prevk;
      if (nh > mMaxRawHashCount) {
        mMaxRawHashCount = nh;
      }
      // make a little histogram
      if (nh <= 10) {
        mMaxHashCountBins[0]++;
      } else if (nh <= 100) {
        mMaxHashCountBins[1]++;
      } else if (nh <= 1000) {
        mMaxHashCountBins[2]++;
      } else if (nh <= 10000) {
        mMaxHashCountBins[3]++;
      } else {
        mMaxHashCountBins[4]++;
      }
      if (mIndexFilterMethod.keepHash(hash, nh)) {
        prevk = k;
        if (nh > mMaxHashCount) {
          mMaxHashCount = (int) nh;
        }
        //System.err.println("hash nh=" + nh);
        mHashCount1 += nh;
        mHashCount2 += nh * nh;
        mNumValues += nh;
        ++mNumHashes;
      }
    }
    mHashCount0 = mNumHashes;
    mHash.trim(mNumValues);
    mValue.trim(mNumValues);

  }

  @Override
  protected final long compressHash(final long hash) {
    return hash;
  }

  @Override
  protected final long decompressHash(final long upper, final long lower) {
    return lower;
  }

  @Override
  public void search(final long hash, final Finder finder) throws IOException {
    if (mState != IndexState.FROZEN) {
      throw new IllegalStateException();
    }
    ++mSearchCount;

    if (!mHashVector.get(hash)) {
      ++mBitVectorMissCount;
      return;
    }
    ++mBitVectorHitCount;

    final long start = position(hash);
    final long compressedHash = compressHash(hash);
    final long low = mInitialPosition.get(start);
    final long high = mInitialPosition.get(start + 1);
    assert low <= high; // : "i=" + i + " low=" + low + " high=" + high;
    final long found = SearchUtils.binarySearch(mHash, low, high - 1, compressedHash);
    if (found < 0) {
      ++mMissCount;
      return;
    }
    ++mInitialHitCount;
    ++mFinalHitCount;
    assert compressedHash == mHash.get(found) : found;
    assert hash == decompressHash(start, mHash.get(found)) : found;
    //System.out.println("found = " + found + " qkey=" + hash + " start=" + start + " skey=" + mHash.get(found));
    long i = found - 1;
    while (i >= low && mHash.get(i) == compressedHash) {
      --i;
    }
    long j = i + 1;
    while (j < high && mHash.get(j) == compressedHash && finder.found(mValue.get(j))) {
      ++mTotalFinds;
      ++j;
    }
  }

  @Override
  public void scan(FinderHashValue finder) throws IOException {
    if (mState != IndexState.FROZEN) {
      throw new IllegalStateException();
    }
    for (long l = 0; l < mNumValues; ++l) {
      final long hash = mHash.get(l);
      final long value = mValue.get(l);
      finder.found(hash, value);
    }
  }

  @Override
  public int count(final long hash) {
    try {
      final CountingFinder countingFinder = new CountingFinder();
      search(hash, countingFinder);
      final long count = countingFinder.getCount();
      assert count <= Integer.MAX_VALUE;
      return (int) count;
    } catch (final IOException e) {
      throw new RuntimeException(e); // CountingFinder has no IOException
    }
  }

  private long find(final long hash) {
    if (mState != IndexState.FROZEN) {
      throw new IllegalStateException();
    }
    ++mSearchCount;

    if (!mHashVector.get(hash)) {
      //assert mInitialPosition.binarySearch(hash) < 0 : "hash=" + hash + LS + toString();
      ++mBitVectorMissCount;
      return -1;
    }
    ++mBitVectorHitCount;
    final long found = binarySearch(hash);
    if (found < 0) {
      ++mMissCount;
    } else {
      ++mInitialHitCount;
    }
    //assert found < 0 || found >= 0 && found < mNumValues;
    return found;
  }

  @Override
  public boolean contains(final long hash) {
    return find(hash) >= 0;
  }

  @Override
  public long first(long hash) {
    long index = find(hash);
    if (index < 0) {
      return index;
    }
    while (true) {
      --index;
      if (index < 0 || getHash(index) != hash) {
        ++index;
        break;
      }
    }
    return index;
  }

  @Override
  public final long getHash(final long found) {
    assert found >= 0 : found;
    return mHash.get(found);
  }

  @Override
  public final long getValue(final long found) {
    return mValue.get(found);
  }

  @Override
  void perfString(final StringBuilder sb) {
    perc(sb, mSearchCount, mSearchCount, "Total search calls");
    perc(sb, mBitVectorMissCount, mSearchCount, "Misses detected by bit vector");
    perc(sb, mBitVectorHitCount, mSearchCount, "Hits detected by bit vector");
    perc(sb, mMissCount, mSearchCount, "Misses on probing hash index");
    perc(sb, mInitialHitCount, mSearchCount, "Hits on probing hash index");
    perc(sb, mFinalHitCount, mSearchCount, "Hits where some entry found");
    sb.append(LS);
    sb.append(LS);
    LONG_FORMAT.format(sb, mTotalFinds).append(" Total finds").append(LS);
    if (mSearchCount != 0) {
      sb.append("    ");
      REAL_FORMAT.format(sb, mTotalFinds / (double) mSearchCount).append(" finds / search").append(LS);
    }
    sb.append(LS);
  }

  @Override
  public boolean integrity() {
    super.integrity();
    Exam.assertFalse(mParams.compressHashes());
    if (mState == IndexState.FROZEN) {
      Exam.assertTrue(mSearchCount >= 0);
      Exam.assertTrue(mBitVectorMissCount >= 0);
      Exam.assertTrue(mBitVectorHitCount >= 0);
      Exam.assertTrue(mMissCount >= 0);
      Exam.assertTrue(mInitialHitCount >= 0);
      Exam.assertTrue(mFinalHitCount >= 0);
      Exam.assertTrue(mInitialHitCount >= 0);
      Exam.assertTrue(mFinalHitCount >= 0);
      Exam.assertTrue(mTotalFinds >= mFinalHitCount);
      Exam.assertEquals(mSearchCount, mBitVectorHitCount + mBitVectorMissCount);
      Exam.assertEquals(mBitVectorHitCount, mMissCount + mInitialHitCount);
      Exam.assertTrue(mTotalFinds >= mFinalHitCount);
      return true;
    } else {
      Exam.assertEquals(0, mSearchCount);
      Exam.assertEquals(0, mBitVectorMissCount);
      Exam.assertEquals(0, mBitVectorHitCount);
      Exam.assertEquals(0, mMissCount);
      Exam.assertEquals(0, mInitialHitCount);
      Exam.assertEquals(0, mFinalHitCount);
      Exam.assertEquals(0, mTotalFinds);
    }
    return true;
  }

  @Override
  public boolean globalIntegrity() {
    integrity();
    for (long l = 0; l < mNumHashes; ++l) {
      Exam.assertTrue(mHashBits == 64 || mHash.get(l) >>> mHashBits == 0);
    }
    if (mState == IndexState.FROZEN) {
      //testing after freeze which includes sorting and construction of
      //overflow, bit vector and initial pointer array
      //hashes sorted
      Exam.assertTrue(ArrayUtils.isSorted(mHash, mNumHashes));
      Exam.assertTrue(ArrayUtils.isSorted(mInitialPosition, 0, mInitialPositionLength));
      for (long l = 0; l < mNumValues; ++l) {
        Exam.assertTrue(mValue.get(l) >= 0);
      }
    } else {
      //testing before freezing, hash array will be unsorted
      //overflow, bitvector and initial pointer are unused
      for (long l = 0; l < mInitialHashes; ++l) {
        Exam.assertTrue(mValue.get(l) >= 0);
      }
    }
    return true;
  }
}
