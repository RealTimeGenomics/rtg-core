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

import com.rtg.util.array.ExtensibleIndex;
import com.rtg.util.array.byteindex.ByteChunks;
import com.rtg.util.array.longindex.LongChunks;

/**
 * Basic bit packing compression.
 */
public class ByteBaseCompression implements ByteCompression {

  private final ByteArray mBytes;
  private final ByteChunks mByteChunks;
  private final ExtensibleIndex mPointers;
  private long mTotalSize;
  private boolean mFrozen;

  /**
   * Basic bit packing for ranges that use less than or equal to 8 bits (no compression on 8 bits).
   * @param range the range of values that can be held.
   */
  //TODO eventually should not use this publicly
  public ByteBaseCompression(int range) {
    assert range <= 256 && range > 0;
    final int minBits = CompressedByteArray.minBits(range);
    if (minBits == 8) {
      mByteChunks = new ByteChunks(0);
      mBytes = null;
    } else {
      mBytes = new CompressedByteArray(0, range, true);
      mByteChunks = null;
    }
    mPointers = new LongChunks(1);
    mPointers.set(0, 0);
    mFrozen = false;
  }

  /**
   * Constructor to hold pre-existing data.
   * This is used by SDF reading.
   * Pointers are 0 based with pointer(i+1) being the exclusive end.
   * @param data the byte array data.
   * @param pointers the pointers into the byte array.
   */
  public ByteBaseCompression(ByteArray data, ExtensibleIndex pointers) {
    mByteChunks = null;
    mBytes = data;
    mPointers = pointers;
    mFrozen = true;
  }

  @Override
  public void add(byte[] buffer, int offset, int length) {
    if (mFrozen) {
      throw new RuntimeException("Adding to a frozen ByteCompression");
    }
    if (mBytes != null) {
      mBytes.set(mTotalSize, buffer, offset, length);
    } else {
      mByteChunks.extendBy(length);
      mByteChunks.copyBytes(buffer, offset, mTotalSize, length);
    }
    mTotalSize += length;
    mPointers.append(mTotalSize);
  }

  /**
   * @param index of a block (0 based).
   * @return the length of the block.
   */
  public int length(long index) {
    return (int) (mPointers.get(index + 1) - mPointers.get(index));
  }


  @Override
  public void get(byte[] buffer, long index, int offset, int length) {
    if (mBytes != null) {
      mBytes.get(buffer, mPointers.get(index) + offset, length);
    } else {
      mByteChunks.getBytes(buffer, 0, mPointers.get(index) + offset, length);
    }
  }

  @Override
  public void freeze() {
    mPointers.trim();
    if (mByteChunks != null) {
      mByteChunks.trim();
    }
    mFrozen = true;
  }

  @Override
  public long bytes() {
    return mPointers.bytes() + (mByteChunks != null ? mByteChunks.bytes() : mBytes.bytes());
  }
}
