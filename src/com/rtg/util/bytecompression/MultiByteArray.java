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
package com.rtg.util.bytecompression;

import java.io.IOException;
import java.io.InputStream;
import java.util.Arrays;

import com.rtg.util.diagnostic.Diagnostic;

/**
 * Distributes values across multiple arrays.
 */
public class MultiByteArray extends ByteArray {
  private final int mBits;
  private final int mChunkSize;
  private final long mMask;
  private byte[][] mData;
  private long mTotalSize;
  private long mCurrentLength;

  /**
   * Constructor
   * @param size the number of values that can be stored
   */
  public MultiByteArray(final long size) {
    this(size, 27);
  }

  /**
   * Constructor
   * @param size the number of values that can be stored
   * @param bits the number of bits of addressing used per array
   */
  public MultiByteArray(final long size, int bits) {
    assert 1 < bits;
    assert bits <= 30;
    mBits = bits;
    mChunkSize = 1 << mBits;
    mMask = mChunkSize - 1;
    mData = new byte[(int) (size >> mBits) + ((size & mMask) == 0 ? 0 : 1)][];
    int i = 0;
    long x = size;
    mTotalSize = 0;
    while (x > 0) {
      final int chunkSize = x > mChunkSize ? mChunkSize : (int) x;
      Diagnostic.developerLog("MultiByteArray allocating " + chunkSize
          + " bytes (block " + (i + 1) + " of " + mData.length + ")");
      mData[i++] = new byte[chunkSize];
      x -= chunkSize;
      mTotalSize += chunkSize;
    }
    mCurrentLength = size;
  }

  /**
   * Ensure there is enough space to store at least given number of entries.
   * Internally there may be room for more.
   * @param newMax new length of array.
   */
  public void extendTo(long newMax) {

    if (newMax < mCurrentLength) {
      throw new IllegalArgumentException("" + newMax + " is less than current length of: " + mCurrentLength);
    }
    while (mTotalSize < newMax) {
      final long i = mTotalSize >>> mBits;
      if (i >= mData.length) {
        final long newSize = (mData.length + 1) * 2;
        if (newSize > Integer.MAX_VALUE) {
          //I wonder what will be the earliest date it will be possible to get to here
          throw new RuntimeException("Attempting to allocate too large a chunk array. newSize=" + newSize);
        }
        mData = Arrays.copyOf(mData, (int) newSize);
      }
      final int ii = (int) i;
      if (mData[ii] == null) {
        mData[ii] = new byte[mChunkSize];
        mTotalSize += mChunkSize;
      } else {
        //short subarray
        final byte[] newSubArray = new byte[mChunkSize];
        final byte[] arr = mData[ii];
        final int lenArr = arr.length;
        System.arraycopy(arr, 0, newSubArray, 0, lenArr);
        mTotalSize += mChunkSize - lenArr;
        mData[ii] = newSubArray;
      }
    }
    mCurrentLength = newMax;
  }

  @Override
  public byte get(long offset) {
    final int block = (int) (offset >> mBits);
    final int blockpos = (int) (offset & mMask);
    return mData[block][blockpos];
  }

  @Override
  public void get(byte[] dest, long offset, int count) {
    int block = (int) (offset >> mBits);
    int blockpos = (int) (offset & mMask);
    int destpos = 0;
    int todo = count;

    while (todo > 0) {
      final int amountToCopy = Math.min(todo, mChunkSize - blockpos);
      System.arraycopy(mData[block], blockpos, dest, destpos, amountToCopy);
      destpos += amountToCopy;
      todo -= amountToCopy;
      block++;
      blockpos = 0;
    }
  }

  @Override
  public void set(long offset, byte value) {
    final int block = (int) (offset >> mBits);
    final int blockpos = (int) (offset & mMask);
    mData[block][blockpos] = value;
  }

  @Override
  public void set(long offset, byte[] buffer, int count) {
    set(offset, buffer, 0, count);
  }

  @Override
  public void set(long offset, byte[] src, int bOffset, int count) {
    int block = (int) (offset >> mBits);
    int blockpos = (int) (offset & mMask);
    int srcpos = bOffset;
    int todo = count;

    while (todo > 0) {
      final int amountToCopy = Math.min(todo, mChunkSize - blockpos);
      System.arraycopy(src, srcpos, mData[block], blockpos, amountToCopy);
      srcpos += amountToCopy;
      todo -= amountToCopy;
      block++;
      blockpos = 0;
    }
  }

  void load(final InputStream stream, final long offset, final int count) throws IOException {
    int read = 0;
    int i = (int) (offset >> mBits);
    int j = (int) (offset & mMask);
    while (read < count) {
      read += stream.read(mData[i], j, Math.min(mChunkSize - j, count - read));
      i++;
      j = 0;
    }
  }

  @Override
  public long length() {
    return mCurrentLength;
  }
}
