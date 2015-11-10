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
import java.util.ConcurrentModificationException;

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.index.params.CreateParams;
import com.rtg.util.IORunnable;
import com.rtg.util.LongUtils;
import com.rtg.util.SimpleThreadPool;
import com.rtg.util.array.ArrayUtils;
import com.rtg.util.array.IndexSorter;
import com.rtg.util.array.Swapper;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.OneShotTimer;
import com.rtg.util.integrity.Exam;

/**
 * Builds and searches a fast index. Two-pass indexes allow for hash compression
 * but require an additional pass that accumulates statistics about the
 * distribution of hash values. Each pass must call add with the same data,
 * followed by a freeze.
 *
 */
@TestClass({"com.rtg.index.IndexCompressedTest", "com.rtg.index.IndexCompressedExtendedTest"})
public class IndexCompressed extends IndexBase implements IndexExtended {

  /** Mask to get the lower bits of a hash value */
  private final long mCompressHashMask;
  private final int mExcessBits;
  private final long mExcessBitsMask;
  private final int mExcessShift;

  /**
   * Constructs an empty index.
   *
   * @param indexParams holds all the values needed for constructing the index.
   * @param threshold maximum repeat frequency threshold - default
   *        <code>Integer.MAX_VALUE</code> if null.
   * @param proportionalThreshold Whether the frequency threshold should be calculated from index data rather than as a parameter.
   * @param maxThreshold when using proportional threshold don't exceed this repeat frequency
   * @param minThreshold when using proportional threshold don't go below this repeat frequency
   * @param numberThreads number of threads appropriate for parallel execution.
   */
  public IndexCompressed(final CreateParams indexParams, final Integer threshold, boolean proportionalThreshold, int maxThreshold, int minThreshold, final int numberThreads) {
    super(indexParams, threshold, proportionalThreshold, maxThreshold, minThreshold, numberThreads);
    mCompressHashMask = LongUtils.longMask(mSRBits);
    mExcessBits = mHashBits > Long.SIZE ? mHashBits - Long.SIZE : 0;
    mExcessBitsMask = LongUtils.longMask(mExcessBits);
    if (mHashBits < 1) {
      throw new RuntimeException("Hash bits set to invalid value bits=" + mHashBits);
    }
    mExcessShift = mInitialPointerBits -  mExcessBits;
    assert mExcessShift >= 0 : "initialPointerBits=" + mInitialPointerBits + " excessBits=" + mExcessBits + " excessShift=" + mExcessShift;
    //System.err.println("hashBits=" + mHashBits + " initialPointerBits=" + mInitialPointerBits + " excessBits=" + mExcessBits + " excessShift=" + mExcessShift);
    //System.err.println("hashMask=" + Utils.toBits(mCompressHashMask, 64));
  }

  @Override
  public final void add(final long hash, final long value) {
    //System.err.println(mState + " add(" + hash + ", " + value + ")");
    if (mState == IndexState.PRE_ADD) {
      final long i = position(hash) + 2;
      mInitialPosition.set(i, mInitialPosition.get(i) + 1);
    } else if (mState == IndexState.FROZEN) {
      throw new IllegalStateException();
    } else {
      assert value >= 0;
      final long upperh = position(hash) + 1;
      final long next = mInitialPosition.get(upperh);
      final long lowerh = compressHash(hash);
      assert next <= mInitialPosition.get(upperh + 1) : "next=" + next + " upperh=" + upperh + " get=" + mInitialPosition.get(upperh + 1);
      mHash.set(next, lowerh);
      mValue.set(next, value);
      mInitialPosition.set(upperh, next + 1);
    }
  }

  @Override
  public void add(long[] hash, long value) {
    if (mState == IndexState.PRE_ADD) {
      final long i = position(hash) + 2;
      mInitialPosition.set(i, mInitialPosition.get(i) + 1);
    } else if (mState == IndexState.FROZEN) {
      throw new IllegalStateException();
    } else {
      assert value >= 0;
      final long upperh = position(hash) + 1;
      final long next = mInitialPosition.get(upperh);
      final long lowerh = compressHash(hash);
      assert next <= mInitialPosition.get(upperh + 1) : "next=" + next + " upperh=" + upperh + " get=" + mInitialPosition.get(upperh + 1);
      mHash.set(next, lowerh);
      mValue.set(next, value);
      mInitialPosition.set(upperh, next + 1);
    }
  }

  @Override
  public void freeze() {
    if (mState == IndexState.PRE_ADD) {
      //System.err.println("freeze1 start");
      //System.err.println(mInitialPosition.toString());
      long sum = 0;
      assert mInitialPosition.get(0) == 0 && mInitialPosition.get(1) == 0;
      for (long i = 1; i < mInitialPositionLength; i++) {
        final long tmp = mInitialPosition.get(i);
        assert tmp >= 0;
        sum += tmp;
        mInitialPosition.set(i, sum);
      }
      mInitialPosition.set(mInitialPositionLength - 1, sum);
      //System.err.println("freeze1 end");
      //System.err.println(mInitialPosition.toString());
      if (sum > mHash.length()) {
        throw new RuntimeException("Too many items pre-added:" + sum + " > " + mHash.length());
      }
      mInitialHashes = sum;
      mState = IndexState.ADD;
    } else if (mState == IndexState.FROZEN) {
      throw new IllegalStateException("Index closed twice");
    } else {
      final OneShotTimer sort = new OneShotTimer("Index_sort");
      //System.err.println("freeze2 start");
      //System.err.println(mInitialPosition.toString());
      final Swapper swapper = new Swapper(mHash, mValue);
      final int numThreads = mNumberThreads;
      final long batchSize = (mInitialPositionLength - 2) / numThreads;
      if ((numThreads == 1) || (batchSize <= 1)) {
        // sort each segment of mHash (and mValue), rather than the whole array
        for (long i = 0; i < mInitialPositionLength - 2; i++) {
          final long lo = mInitialPosition.get(i);
          final long hi = mInitialPosition.get(i + 1);
          assert lo <= hi : "lo=" + lo + " hi=" + hi + " i=" + i;
          IndexSorter.sort(mHash, lo, hi - lo, swapper);
        }
      } else {
        if (mHash.safeFromWordTearing()) {
          final SimpleThreadPool sp = new SimpleThreadPool(numThreads, "Building", true);
          for (long i = 0; i < numThreads; i++) {
            sp.execute(new SwapThread(batchSize * i, (i == numThreads - 1) ? (mInitialPositionLength - 2) : (batchSize * (i + 1)), swapper));
          }
          try {
            sp.terminate();
          } catch (final IOException e) {
            throw new IllegalStateException("Index sorting should not throw IOException", e);
          }
        } else {
          throw new ConcurrentModificationException("Index implementation is not safe from concurrent update of adjacent values");
        }

      }
      //System.err.println("freeze2 end");
      //System.err.println(mInitialPosition.toString());
      sort.stopLog();
      //fill overflow table
      final OneShotTimer over = new OneShotTimer("Index_pointer");
      createPointers();
      over.stopLog();

      if (mHashVector != null) {
        final OneShotTimer bitv = new OneShotTimer("Index_bitVector");
        //set bit vector
        // we loop through mInitialPosition to get upper bits of hash
        long low = 0;
        for (long upper = 0; upper < mInitialPositionLength - 2; upper++) {
          final long hi = mInitialPosition.get(upper + 1);
          assert low <= hi;
          if (low < hi) {
            long lastHash = mHash.get(low);
            final long deh0 = mExcessBitsMask == 0 ? decompressHash(upper, lastHash) : decompressHashExtended(upper, lastHash)[0];
            mHashVector.set(deh0);
            for (long i = low + 1; i < hi; i++) {
              final long lower = mHash.get(i);
              if (lower == lastHash) {
                continue;
              }
              final long deh = mExcessBitsMask == 0 ? decompressHash(upper, lower) : decompressHashExtended(upper, lower)[0];
              mHashVector.set(deh);
              lastHash = lower;
            }
          }
          low = hi;
        }
        assert low == mNumValues; // we should have processed each hash once
        bitv.stopLog();
      }

      checkRepeatFrequency();

      mState = IndexState.FROZEN;
      assert globalIntegrity();
    }
  }

  private class SwapThread implements IORunnable {

    final long mStart;
    final long mEnd;
    final Swapper mSwapper;

    SwapThread(final long start, final long end, final Swapper swapper) {
      mStart = start;
      mEnd = end;
      mSwapper = swapper;
    }

    @Override
    public void run() {
      //System.err.println("mt sort " + mStart + " to " + mEnd);
      for (long i = mStart; i < mEnd; i++) {
        final long lo = mInitialPosition.get(i);
        final long hi = mInitialPosition.get(i + 1);
        assert lo <= hi : "mStart=" + mStart + " mEnd=" + mEnd + " lo=" + lo + " hi=" + hi + " i=" + i;
        IndexSorter.sort(mHash, lo, hi - lo, mSwapper);
      }

    }
  }

  private int determineCutoffFrequency() {
    final OneShotTimer over = new OneShotTimer("Index_Frequency");
    int[] freqDist = new int[1024];
    SparseFrequencyHistogram freqHist = new SparseFrequencyHistogram();
    long mergeCount = 0;
    int numUsed = 0;
    long lo = 0;
    for (long p = 0; p < mInitialPositionLength - 2; p++) {
      final long hi = mInitialPosition.get(p + 1);
      for (long i = lo; i < hi;) {
        final long hash = mHash.get(i);
        int freq = 1;
        i++;
        while (i < hi && hash == mHash.get(i)) {
          i++;
          freq++;
        }
        if (numUsed >= freqDist.length) {
          if (freqDist.length == MAX_FREQ_DIST_SIZE) {
            freqHist = SparseFrequencyHistogram.merge(freqHist, SparseFrequencyHistogram.fromIndividualFrequencies(freqDist, numUsed));
            mergeCount++;
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
      lo = hi;
    }
    freqHist = SparseFrequencyHistogram.merge(freqHist, SparseFrequencyHistogram.fromIndividualFrequencies(freqDist, numUsed));
    mergeCount++;
    Diagnostic.developerLog("IndexCompressed " + mergeCount + " Frequency Histogram Merges");
    final int ret = determineCutoffFrequency(freqHist, mThreshold);
    over.stopLog();
    return ret;
  }

  @Override
  protected void compact() {
    mMaxHashCount = 0;
    final long threshold;
    if (mUseProportionalThreshold) {
      threshold = determineCutoffFrequency();
    } else {
      threshold = mThreshold;
    }
    //System.err.println("compact");
    long prevk = 0;
    long lo = 0;
    for (long p = 0; p < mInitialPositionLength - 2; p++) {
      final long hi = mInitialPosition.get(p + 1);
      for (long i = lo; i < hi;) {
        long k = prevk;
        final long hash = mHash.get(i);
        //optimize
        while (i < hi && hash == mHash.get(i)) {
          mHash.set(k, hash);
          mValue.set(k, mValue.get(i));
          k++;
          i++;
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

        if (nh <= threshold) {
          prevk = k;
          //System.err.println("hash nh=" + nh);

          if (nh > mMaxHashCount) {
            mMaxHashCount = (int) nh;
          }
          mHashCount1 += nh;
          mHashCount2 += nh * nh;
          mNumValues += nh;
          mNumHashes++;
        }
      }
      lo = hi;
      mInitialPosition.set(p + 1, prevk);
      final long bucketCnt = prevk - mInitialPosition.get(p);
      if (bucketCnt > 0) {
        mBucketCount0++;
        mBucketCount1 += bucketCnt;
        mBucketCount2 += bucketCnt * bucketCnt;
      }
    }
    mHashCount0 = mNumHashes;
    assert mHashCount1 == mNumValues; // : mHashCount1 + ":" + mInitialHashes;
    mInitialPosition.set(mInitialPositionLength - 1, mInitialPosition.get(mInitialPositionLength - 2));
    mHash.trim(mNumValues);
    mValue.trim(mNumValues);
  }

  /**
   * Extract the lower bits from a hash. That is, the bits that are NOT used as
   * an index into <code>mInitialPosition</code>.
   *
   * @param hash being sought.
   * @return lower bits.
   */
  @Override
  protected final long compressHash(final long hash) {
    assert mExcessBits == 0;
    return hash & mCompressHashMask;
  }

  /**
   * Extract the bits used for an extended hash. That is, the bits that are NOT used as
   * an index into <code>mInitialPosition</code>.
   *
   * @param hash being sought.
   * @return used bits.
   */
  protected final long compressHash(final long[] hash) {
    return hash[0] & mCompressHashMask;
  }

  /**
   * Reconstruct a hash value from two components.
   *
   * @param upper an index into <code>mInitialPosition</code>.
   * @param lower the other bits.
   * @return a full hash
   */
  @Override
  protected final long decompressHash(final long upper, final long lower) {
    assert 0 <= lower && lower <= mCompressHashMask;
    assert 0 <= upper && upper < mInitialPositionLength - 1;
    return (upper - mIncrement) << mSRBits | lower;
  }

  /**
   * Reconstruct a hash value from two components.
   *
   * @param upper an index into <code>mInitialPosition</code>.
   * @param lower the other bits.
   * @return a full hash
   */
  protected final long[] decompressHashExtended(final long upper, final long lower) {
    //assert 0 <= upper && upper < mInitialPositionLength - 1;
    //System.err.println("sRbits=" + mSRBits + " excessBits=" + mExcessBits);
    //System.err.println(" upper=" + Utils.toBits(upper, 64));
    //System.err.println(" lower=" + Utils.toBits(lower, 64));
    final long up = upper - mIncrement;
    if (mHashBits > Long.SIZE) {
      final long[] dec = new long[2];
      dec[1] = up >>> mExcessShift;
    if (mExcessShift == 0) {
      dec[0] = lower;
    } else {
      assert mExcessShift > 0;
      dec[0] = (up << (Long.SIZE - mExcessShift)) | lower;
    }
    return dec;
    }
    assert mExcessBits == 0;
    final long[] dec = new long[1];
    dec[0] = up << mSRBits | lower;
    return dec;
  }

  @Override
  protected  long position(final long hash) {
    final long k = (hash >> mSRBits) + mIncrement;
    assert k >= 0 && k < mInitialPositionLength - 1 : "k=" + k + " hash=" + hash + " srbits=" + mSRBits;
    return k;
  }

  /**
   * Convert from a hash to an index in <code>mInitialPosition</code>.
   * @param hash being sought.
   * @return index into <code>mInitialPosition</code>.
   */
  protected final long position(final long[] hash) {
    if (mHashBits <= Long.SIZE) {
      return position(hash[0]);
    }
    final long h1 = hash[hash.length - 1];
    if (mExcessShift == 0) {
      return h1;
    }
    final long h2 = hash[hash.length - 2];
    return (h1 << mExcessShift) | h2 >>> (Long.SIZE - mExcessShift);
  }


  @Override
  public void search(final long hash, final Finder finder) throws IOException, IllegalStateException {
    if (mState != IndexState.FROZEN) {
      throw new IllegalStateException();
    }

    if (!mHashVector.get(hash)) {
      return;
    }

    final long start = position(hash);
    final long compressedHash = compressHash(hash);
    final long low = mInitialPosition.get(start);
    final long high = mInitialPosition.get(start + 1);
    assert low <= high; // : "i=" + i + " low=" + low + " high=" + high;
    final long found = SearchUtils.binarySearch(mHash, low, high - 1, compressedHash);
    if (found < 0) {
      return;
    }
    assert compressedHash == mHash.get(found) : found;
    assert hash == decompressHash(start, mHash.get(found)) : found;
    //System.out.println("found = " + found + " qkey=" + hash + " start=" + start + " skey=" + mHash.get(found));
    long i = found - 1;
    while (i >= low && mHash.get(i) == compressedHash) {
      i--;
    }
    long j = i + 1;
    while (j < high && mHash.get(j) == compressedHash && finder.found(mValue.get(j))) {
      j++;
    }
  }

  @Override
  public void search(long[] hash, Finder finder) throws IOException, IllegalStateException {

    if (mState != IndexState.FROZEN) {
      throw new IllegalStateException();
    }

    if (!mHashVector.get(hash[0])) {
      return;
    }

    final long start = position(hash);
    final long compressedHash = compressHash(hash);
    final long low = mInitialPosition.get(start);
    final long high = mInitialPosition.get(start + 1);
    assert low <= high; // : "i=" + i + " low=" + low + " high=" + high;
    final long found = SearchUtils.binarySearch(mHash, low, high - 1, compressedHash);
    if (found < 0) {
      return;
    }
    assert compressedHash == mHash.get(found) : found;
    assert Arrays.equals(hash, decompressHashExtended(start, mHash.get(found))) : found;
    //System.out.println("found = " + found + " qkey=" + hash + " start=" + start + " skey=" + mHash.get(found));
    long i = found - 1;
    while (i >= low && mHash.get(i) == compressedHash) {
      i--;
    }
    long j = i + 1;
    while (j < high && mHash.get(j) == compressedHash && finder.found(mValue.get(j))) {
      j++;
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
      throw new RuntimeException("Unpossible: " + e.getMessage()); // CountingFinder has no IOException
    }
  }

  @Override
  public int count(long[] hash) {
    throw new UnsupportedOperationException();
  }

  @Override
  public void scan(FinderHashValue finder) throws IOException, IllegalStateException {
    if (mState != IndexState.FROZEN) {
      throw new IllegalStateException();
    }
    long lo = 0;
    for (long p = 0; p < mInitialPositionLength - 2; p++) {
      final long hi = mInitialPosition.get(p + 1);
      for (long i = lo; i < hi; i++) {
        final long hash0 = mHash.get(i);
        final long hash = decompressHash(p, hash0);
        final long value = mValue.get(i);
        finder.found(hash, value);
      }
      lo = hi;
    }
  }

  @Override
  public void scanAll(FinderHashValueExtended finder) throws IOException, IllegalStateException {
    if (mState != IndexState.FROZEN) {
      throw new IllegalStateException();
    }
    long lo = 0;
    for (long p = 0; p < mInitialPositionLength - 2; p++) {
      final long hi = mInitialPosition.get(p + 1);
      for (long i = lo; i < hi; i++) {
        final long hash0 = mHash.get(i);
        final long[] hash = decompressHashExtended(p, hash0);
        final long value = mValue.get(i);
        finder.found(hash, value);
      }
      lo = hi;
    }
  }


  private long find(final long hash) {
    if (mState != IndexState.FROZEN) {
      throw new IllegalStateException();
    }
    if (!mHashVector.get(hash)) {
      //assert mInitialPosition.binarySearch(hash) < 0 : "hash=" + hash + LS + toString();
      return -1;
    }
    final long start = position(hash);
    final long low = mInitialPosition.get(start);
    final long high = mInitialPosition.get(start + 1);
    assert low <= high; // : "i=" + i + " low=" + low + " high=" + high;
    final long compressedHash = compressHash(hash);
    return SearchUtils.binarySearch(mHash, low, high - 1, compressedHash);
  }

  @Override
  public boolean contains(final long hash) {
    return find(hash) >= 0;
  }

  @Override
  public long first(long hash) throws IllegalStateException {
    long index = find(hash);
    if (index < 0) {
      return index;
    }
    while (true) {
      index--;
      if (index < 0 || getHash(index) != hash) {
        index++;
        break;
      }
    }
    return index;
  }

  private long find(long[] hash) {
    if (mState != IndexState.FROZEN) {
      throw new IllegalStateException();
    }
    if (!mHashVector.get(hash[0])) {
      //assert mInitialPosition.binarySearch(hash) < 0 : "hash=" + hash + LS + toString();
      return -1;
    }
    final long start = position(hash);
    final long low = mInitialPosition.get(start);
    final long high = mInitialPosition.get(start + 1);
    assert low <= high; // : "i=" + i + " low=" + low + " high=" + high;
    final long compressedHash = compressHash(hash);
    return SearchUtils.binarySearch(mHash, low, high - 1, compressedHash);
  }

  @Override
  public boolean contains(final long[] hash) {
    return find(hash) >= 0;
  }

  @Override
  public long first(long[] hash) throws IllegalStateException {
    long index = find(hash);
    if (index < 0) {
      return index;
    }
    while (true) {
      index--;
      if (index < 0 || !Arrays.equals(getHashExtended(index), hash)) {
        index++;
        break;
      }
    }
    return index;
  }

  @Override
  public final long getHash(final long found) {
    final long index = SearchUtils.bracketSearch(mInitialPosition, 0, mInitialPositionLength - 1, found);
    assert index >= 0 : "found=" + found + " index=" + index;
    return decompressHash(index, mHash.get(found));
  }

  @Override
  public long[] getHashExtended(long found) {
    final long index = SearchUtils.bracketSearch(mInitialPosition, 0, mInitialPositionLength - 1, found);
    assert index >= 0 : "found=" + found + " index=" + index;
    return decompressHashExtended(index, mHash.get(found));
  }

  @Override
  public long getValue(long found) {
    return mValue.get(found);
  }

  @Override
  void perfString(final StringBuilder sb) {
    sb.append("Performance statistics not available.").append(LS);
  }

  @Override
  public boolean integrity() {
    super.integrity();
    if (mState == IndexState.FROZEN) {
      Exam.assertEquals(mInitialPositionLength + ":" + mInitialPosition.toString(), mInitialPosition.get(mInitialPositionLength - 1), mInitialPosition.get(mInitialPositionLength - 2));
    }
    Exam.assertTrue(mExcessBits == 0 || mExcessBits == mHashBits - Long.SIZE);
    return true;
  }

  @Override
  public boolean globalIntegrity() {
    integrity();
    if (mHashBits < Long.SIZE) {
      for (long l = 0; l < mNumHashes; l++) {
        Exam.assertTrue(mHash.get(l) >>> mHashBits == 0);
      }
    }
    if (mState == IndexState.FROZEN) {
      //testing after freeze which includes sorting and construction of
      //overflow, bit vector and initial pointer array
      for (long i = 1; i < mInitialPositionLength - 1; i++) {
        final long lo = mInitialPosition.get(i - 1);
        final long hi = mInitialPosition.get(i);
        Exam.assertTrue(ArrayUtils.isSorted(mHash, lo, hi));
      }
      Exam.assertTrue(ArrayUtils.isSorted(mInitialPosition, 0, mInitialPositionLength));
      //Initial position
      final long first = mInitialPosition.get(0);
      Exam.assertTrue("first=" + first, first >= 0 && first <= mNumHashes);
      for (long i = 1; i < mInitialPositionLength; i++) {
        final long lo = mInitialPosition.get(i - 1);
        final long hi = mInitialPosition.get(i);
        Exam.assertTrue(lo <= hi);
        Exam.assertTrue(i + ":" + hi + ":" + mNumValues, 0 <= hi && hi <= mNumValues);
      }
    }
    return true;
  }
}
