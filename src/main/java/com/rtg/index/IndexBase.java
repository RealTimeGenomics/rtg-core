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
package com.rtg.index;

import java.io.PrintStream;

import com.rtg.index.params.CreateParams;
import com.rtg.util.StringUtils;
import com.rtg.util.array.CommonIndex;
import com.rtg.util.array.ExtensibleIndex;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.OneShotTimer;
import com.rtg.util.format.FormatInteger;
import com.rtg.util.format.FormatReal;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;

/**
 * Builds and searches a fast index.
 *
 */
public abstract class IndexBase extends IntegralAbstract implements Index {

  protected static final FormatReal PERC_FORMAT = new FormatReal(5, 1);

  static final int MAX_FREQ_DIST_SIZE = 30000000; //Approximately 100 MB worth by default

  private static final int FORMAT_LENGTH = 30;
  static final FormatInteger LONG_FORMAT = new FormatInteger(FORMAT_LENGTH, true);
  static final FormatReal REAL_FORMAT = new FormatReal(FORMAT_LENGTH - 4, 3);

  protected final CreateParams mParams;

  //The following all deal with the initial pointer data
  protected final CommonIndex mInitialPosition;

  /** Number of bits used for table - it is 2^bits + 1 long. */
  protected final int mInitialPointerBits;

  /** Number of bits to shift the key right in order to extract the index. */
  protected final int mSRBits;

  /** Length of array. */
  protected final long mInitialPositionLength;

  /** Stores the hash of every window in the database */
  protected final ExtensibleIndex mHash;

  /** Values associated with each hash. */
  protected final ExtensibleIndex mValue;

  protected final HashBitVector mHashVector;

  /** Number of bits actually used in the hash code. */
  protected final int mHashBits;

  /**
   * The initial number of hashes before constructing pointers and removing
   * duplicates.
   */
  protected long mInitialHashes = 0;

  /**
   * The number of distinct hashes that have been stored after removing
   * duplicates and setting pointer array.
   */
  protected long mNumHashes = 0;

  /** The number of distinct values <code>mHash</code>. */
  protected long mNumValues = 0;

  /** The maximum hash count */
  protected long mMaxRawHashCount = 0;

  /** Maximum hash count post compact */
  protected int mMaxHashCount = -1;

  /** The maximum hash count in a few bins, &lt; 10, 11-100, 101-1000, 10001-10000, &gt;100000 */
  protected final long[] mMaxHashCountBins = new long[5];

  protected long mHashCount0 = 0;
  protected long mHashCount1 = 0;
  protected long mHashCount2 = 0;

  protected long mBucketCount0 = 0;
  protected long mBucketCount1 = 0;
  protected long mBucketCount2 = 0;

  /**
   * The state the current index is in.
   *
   */
  protected enum IndexState {
    /**
     * Pre-add encompasses both the pre-add and pre-freeze statistic calculation
     * steps
     */
    PRE_ADD,
    /** Hashes are being added */
    ADD,
    /**
     * The system has been frozen - that is hash is sorted, and overflow bit
     * vector and initial pointer have been frozen.
     */
    FROZEN
  }

  protected IndexState mState;

  /**
   * An offset so that negative longs get set as positive positions in the
   * correct order. Only an issue for 64 bit hashes (all others have sign bit ==
   * 0)
   */
  protected final long mIncrement;

  protected final int mNumberThreads;

  protected final IndexFilterMethod mIndexFilterMethod;

  /**
   * Constructs an empty index.
   *
   * @param indexParams holds all the values needed for constructing the index.
   * @param filter the handle to be used for filtering hashes
   * @param numberThreads number of threads appropriate for parallel execution.
   */
  public IndexBase(final CreateParams indexParams, IndexFilterMethod filter, final int numberThreads) {
    assert indexParams.integrity();
    mParams = indexParams;
    mIndexFilterMethod = filter;

    final long size = indexParams.size();
    final int hashBits = indexParams.hashBits();
    Diagnostic.developerLog("hashBits=" + hashBits);
    if (size < 0) {
      throw new RuntimeException("Size must be positive:" + size);
    }
    final OneShotTimer init = new OneShotTimer("Index_initialization");

    mState = mParams.compressHashes() ? IndexState.PRE_ADD : IndexState.ADD;

    mHashBits = mParams.hashBits();
    mNumberThreads = numberThreads;
    //carefully allocate memory
    try {
      mInitialPointerBits = mParams.initialPointerBits();
      mInitialPosition = mParams.initialPosition().createUnsigned();
      mInitialPositionLength = mInitialPosition.length();
      mSRBits = mInitialPointerBits >= mHashBits ? 0 : mHashBits - mInitialPointerBits;
      //System.err.println("initialPointerBits=" + mInitialPointerBits + " hashBits=" + mHashBits + " SRBits=" + mSRBits);

      mHash = mParams.hash().createUnsigned();
      final HashBitHandle bitVector = mParams.bitVector();
      if (bitVector != null) {
        mHashVector = bitVector.create();
      } else {
        mHashVector = null;
      }
      mValue = mParams.value().createUnsigned();

      assert mSRBits >= 0;
      if (mHashBits == 64) {
        mIncrement = -(Long.MIN_VALUE >> mSRBits);
      } else {
        mIncrement = 0;
      }

      init.stopLog();
    } catch (final OutOfMemoryError e) {
      Diagnostic.oomMessage();
      Diagnostic.userLog(e);
      throw e;
    }
  }

  @Override
  public abstract void add(final long hash, final long value);

  private static volatile boolean sHavePrintedRepeatFrequencyWarning = false;
  private static final double REPEAT_FREQUENCY_REPORTING_FRACTION = 0.5;

  private static synchronized void printWarning(int percent, long numValues, double initialHashes) {
    if (!sHavePrintedRepeatFrequencyWarning && numValues < initialHashes * REPEAT_FREQUENCY_REPORTING_FRACTION) {
      sHavePrintedRepeatFrequencyWarning = true;
      Diagnostic.warning(percent + "% of window values were discarded. Increase -w or --repeat-freq to reduce this.");
    }
  }

  protected void checkRepeatFrequency() {
    if (mInitialHashes != 0) {
      final long discarded = mInitialHashes - mNumValues;
      final int percent = (int) (100.0 * discarded / mInitialHashes + 0.5);
      Diagnostic.userLog("Applied index hash filtering (initial " + mInitialHashes + ", now " + mNumValues + ", " + percent + "% of hashes discarded)");
      if (!sHavePrintedRepeatFrequencyWarning) {
        printWarning(percent, mNumValues, mInitialHashes);
      }
    }
  }


  /**
   * Returns the maximum hash count
   *
   * @return the maximum hash count
   */
  public final long maxRawHashCount() {
    return mMaxRawHashCount;
  }

  /**
   * Returns the hash count array
   *
   * @return hash count array
   */
  public String hashCountArray() {
    return "(0-10/100/1000/10000/10001+) " + mMaxHashCountBins[0] + ", "
        + mMaxHashCountBins[1] + ", " + mMaxHashCountBins[2] + ", " + mMaxHashCountBins[3] + ", " + mMaxHashCountBins[4];
  }

  /**
   * Create pointers if actually doing overflow. Throw away hashes that exceed
   * threshold.
   */
  protected final void createPointers() {
    if (mInitialHashes == 0) {
      //special case - easier to return here than get later code correct
      mNumHashes = 0;
      mMaxHashCount = 0;
      return;
    }
    compact();
    Diagnostic.userLog("Hash count histogram: " + hashCountArray());
    Diagnostic.userLog("Maximum hash count: " + maxRawHashCount());
    Diagnostic.userLog("Maximum hash count post compact: " + maxHashCount());
    printStatistics("post compact statistics:");
  }

  private void printStatistics(String head) {
    Diagnostic.developerLog(head);
    Diagnostic.developerLog("Num hashes = " + StringUtils.commas(mNumHashes));
    Diagnostic.developerLog("Num values = " + StringUtils.commas(mNumValues));
    Diagnostic.developerLog("Num initial position length = " + StringUtils.commas(mInitialPositionLength));
    Diagnostic.developerLog("Initial position array size : " + StringUtils.commas(mInitialPosition.bytes()) + " bytes, length = " + StringUtils.commas(mInitialPosition.length()));
    Diagnostic.developerLog("Hash array size : " + StringUtils.commas(mHash.bytes()) + " bytes, length = " + StringUtils.commas(mHash.length()));
    Diagnostic.developerLog("Position array size : " + StringUtils.commas(mValue.bytes()) + " bytes, length = " + StringUtils.commas(mValue.length()));
    if (mHashVector != null) {
      Diagnostic.developerLog("bit vector : " + StringUtils.commas(mHashVector.bytes()) + " bytes, length = " + StringUtils.commas(mHashVector.length()));
    }
  }

  /**
   * Throw away hashes that exceed threshold, using compressed
   * <code>mHash</code>. This also updates the pointers in
   * <code>mInitialPosition</code>.
   */
  protected abstract void compact();

  /**
   * Perform a binary search in the hash array guided by the initial position
   * array.
   *
   * @param hash being sought.
   * @return its position in the hash table, -1 if it is not found.
   */
  protected final long binarySearch(final long hash) {
    final long i = position(hash);
    final long low = mInitialPosition.get(i);
    final long high = mInitialPosition.get(i + 1) - 1;
    assert low <= high + 1; // : "i=" + i + " low=" + low + " high=" + high;
    return SearchUtils.binarySearch(mHash, low, high, hash);
  }

  /**
   * Extract the lower bits from a hash. That is, the bits that are NOT used as
   * an index into <code>mInitialPosition</code>.
   *
   * @param hash being sought.
   * @return lower bits.
   */
  protected abstract long compressHash(final long hash);

  /**
   * Reconstruct a hash value from two components.
   *
   * @param upper an index into <code>mInitialPosition</code>.
   * @param lower the other bits.
   * @return a full hash
   */
  protected abstract long decompressHash(final long upper, final long lower);

  /**
   * Convert from a hash to an index in <code>mInitialPosition</code>. Isn't
   * inlined so we can do the careful assertions.
   *
   * @param hash being sought.
   * @return index into <code>mInitialPosition</code>.
   */
  protected abstract long position(final long hash);

  @Override
  public long numberEntries() {
    if (mState != IndexState.FROZEN) {
      throw new IllegalStateException();
    }
    return mNumValues;
  }

  @Override
  public long numberHashes() {
    if (mState != IndexState.FROZEN) {
      throw new IllegalStateException();
    }
    return mNumHashes;
  }

  @Override
  public long getInitialHashes() {
    return mInitialHashes;
  }

  @Override
  public String infoString() {
    final StringBuilder sb = new StringBuilder();
    infoString(sb);
    return sb.toString();
  }

  @Override
  public long bytes() {
    long totalBytes = 0;
    totalBytes += mHash.bytes();
    totalBytes += mValue.bytes();
    totalBytes += mInitialPosition.bytes();
    if (mHashVector != null) {
      totalBytes += mHashVector.bytes();
    }
    return totalBytes;
  }

  void infoString(final StringBuilder sb) {
    sb.append("Memory Usage\tbytes\tlength").append(LS);
    long totalBytes = 0;
    sb.append("\t\t").append(StringUtils.commas(mHash.bytes())).append("\t").append(StringUtils.commas(mHash.length())).append("\tHash").append(LS);
    totalBytes += mHash.bytes();

    sb.append("\t\t").append(StringUtils.commas(mValue.bytes())).append("\t").append(StringUtils.commas(mValue.length())).append("\tValue").append(LS);
    totalBytes += mValue.bytes();

    sb.append("\t\t").append(StringUtils.commas(mInitialPosition.bytes())).append("\t").append(StringUtils.commas(mInitialPosition.length())).append("\tInitial Position")
    .append(LS);
    totalBytes += mInitialPosition.bytes();

    if (mHashVector != null) {
      sb.append("\t\t").append(StringUtils.commas(mHashVector.bytes())).append("\t").append(StringUtils.commas(mHashVector.length())).append("\tBit vector").append(LS);
      totalBytes += mHashVector.bytes();
    }

    sb.append("\t\t").append(StringUtils.commas(totalBytes)).append("\t\tTotal bytes").append(LS);
    if (totalBytes != bytes()) {
      throw new RuntimeException("Inconsistency in memory calculation bytes()=" + bytes() + " total=" + totalBytes);
    }

    sb.append(LS);
    sb.append("Hash counts\t0\t1\t2").append(LS);
    sb.append("\t\t").append(mHashCount0).append("\t").append(mHashCount1).append("\t").append(mHashCount2).append(LS);
    sb.append("Bucket counts\t0\t1\t2").append(LS);
    sb.append("\t\t").append(mBucketCount0).append("\t").append(mBucketCount1).append("\t").append(mBucketCount2).append(LS);
  }

  @Override
  public String perfString() {
    final StringBuilder sb = new StringBuilder();
    perfString(sb);
    return sb.toString();
  }

  /**
   * Used to output numbers that are part of hit/miss counts.
   * @param sb buffer where output being placed.
   * @param count number being displayed.
   * @param total total number of which count is a percentage
   * @param msg accompanying message.
   */
  protected static void perc(final StringBuilder sb, final long count, long total, final String msg) {
    LONG_FORMAT.format(sb, count);
    sb.append("  ");
    PERC_FORMAT.format(sb, total == 0 ? 0 : (100.0 * (count / (double) total)));
    sb.append("% ");
    sb.append(msg);
    sb.append(LS);
  }

  abstract void perfString(final StringBuilder sb);

  @Override
  public boolean integrity() {
    Exam.assertTrue(mInitialHashes >= 0);
    Exam.assertTrue(mValue != null);
    if (mHash == null) {
      Exam.assertTrue(false);
    } else {
      Exam.assertEquals(mHash.length(), mValue.length());
    }
    Exam.assertTrue(mInitialPosition != null);
    if (mState == IndexState.FROZEN) {
      Exam.assertTrue(mNumHashes <= mNumValues);
      //initial position
      if (mInitialPosition == null) {
        Exam.assertTrue(false);
      } else {
        Exam.assertEquals(mInitialPositionLength, mInitialPosition.length());
      }
      Exam.assertTrue(mInitialPositionLength == (1 << mInitialPointerBits) + 1 || mInitialPositionLength == (1 << mInitialPointerBits) + 2);
      Exam.assertTrue(mInitialPositionLength > 1);
      Exam.assertTrue(mInitialPointerBits >= mHashBits || mHashBits == mInitialPointerBits + mSRBits);
      Exam.assertTrue(mInitialPointerBits >= 0);
      Exam.assertTrue(mSRBits >= 0);
      Exam.assertTrue(mHashCount1 <= mBucketCount1); //bucket counts are made before compaction
      Exam.assertTrue(0 <= mHashCount0 && mHashCount0 <= mHashCount1 && mHashCount1 <= mHashCount2);
      Exam.assertTrue(mBucketCount0 + ":" + mBucketCount1 + ":" + mBucketCount2, 0 <= mBucketCount0 && mBucketCount0 <= mBucketCount1 && mBucketCount1 <= mBucketCount2);
      return true;
    } else {
      //testing before freezing, hash array will be unsorted
      //overflow, bitvector and initial poonter are unused
      //System.err.println("initial=" + mInitialHashes + " num=" + mNumHashes);
      Exam.assertEquals(0, mNumHashes);
    }
    if (mState == IndexState.FROZEN) {
      Exam.assertTrue(mNumHashes + ":" + mInitialHashes, mNumHashes <= mInitialHashes);
      Exam.assertTrue(mNumValues <= mInitialHashes);
    } else {
      Exam.assertTrue(0 <= mInitialHashes);
    }
    return true;
  }

  @Override
  public final void toString(final StringBuilder sb) {
    sb.append("IndexImplementation").append(LS);
    sb.append("Hash Bits=").append(mHashBits).append(LS);
    sb.append("Initial Number of hashes=").append(mInitialHashes).append(LS);
    sb.append("Number of hash codes=").append(mNumHashes).append(LS);
    infoString(sb);
    final long len = mState == IndexState.FROZEN ? mNumValues : mInitialHashes;
    if (len > 0) {
      sb.append(LS);
      sb.append("\tHash" + "\tValue").append(LS);
      for (long i = 0; i < len; ++i) {
        sb.append("[").append(i).append("]\t");
        sb.append(mHash.get(i)).append("\t");
        sb.append(mValue.get(i)).append(LS);
      }
      sb.append(LS);
    }
    if (mState == IndexState.FROZEN) {
      sb.append("Frozen").append(LS);
      perfString(sb);
      sb.append("Initial position").append(LS);
      sb.append("pointer bits=").append(mInitialPointerBits).append(LS).append("shift bits=").append(mSRBits).append(LS).append("length=").append(mInitialPositionLength).append(LS);
      sb.append(mInitialPosition);
      sb.append(LS);

      if (mHashVector != null) {
        sb.append(LS).append("Bit vector").append(LS);
        mHashVector.toString(sb);
      }
    } else if (mState == IndexState.PRE_ADD) {
      sb.append("Gaseous").append(LS).append(LS);
    } else {
      sb.append("Melted").append(LS).append(LS);
    }
  }

  /**
   * Will write out hashes and values in a format identical for both overflow
   * and no overflow cases. Needed to investigate a bug.
   *
   * @param out where to put output.
   */
  @Override
  public void dumpValues(final PrintStream out) {
    if (mState != IndexState.FROZEN) {
      throw new RuntimeException();
    }
    out.println("Index InitialPosition");
    for (long i = 0; i < mInitialPosition.length(); ++i) {
      out.println("[" + i + "]" + "  " + mInitialPosition.get(i));
    }
    out.println("Index Hash  Values");
    for (long i = 0; i < mNumValues; ++i) {
      out.println("[" + i + "]" + "  " + mHash.get(i) + "  " + mValue.get(i));
    }
  }

  @Override
  public int maxHashCount() {
    if (mMaxHashCount == -1) {
      throw new IllegalStateException("Count only available after compact has been called.");
    }
    return mMaxHashCount;
  }

}
