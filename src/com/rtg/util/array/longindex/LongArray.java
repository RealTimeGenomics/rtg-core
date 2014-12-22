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
package com.rtg.util.array.longindex;

import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;

import com.rtg.util.array.IndexType;

/**
 * Index is implemented using a straight forward long array.
 * This is so that short instances of LongIndex can be as efficient as possible.
 *
 */
public final class LongArray extends LongIndex {

  private final long[] mArray;

  /**
   * This should be called directly only in tests.
   *
   * @param length number of items to be stored.
   */
  public LongArray(final long length) {
    super(length);
    //assert length <= MAX_LENGTH;
    mArray = new long[(int) length];
  }

  private LongArray(long[] data, long length) {
    super(length);
    mArray = data;
  }

  @Override
  public long get(final long index) {
    final int ii = (int) index;
    if (ii != index) {
      throw new IndexOutOfBoundsException("ii=" + ii + " mArrays.length=" + mArray.length + " index=" + index);
    }
    return mArray[ii];
  }

  @Override
  public void set(final long index, final long value) {
    final int ii = (int) index;
    if (ii != index) {
      throw new IndexOutOfBoundsException(String.valueOf(index));
    }
    mArray[ii] = value;
  }

  @Override
  public void swap(long index1, long index2) {
    final int ii1 = (int) index1;
    final int ii2 = (int) index2;
    final long tmp = mArray[ii1];
    mArray[ii1] = mArray[ii2];
    mArray[ii2] = tmp;
  }

  @Override
  public boolean integrity() {
    super.integrity();
    assert mArray == null || mArray.length == mLength;
    return true;
  }

  @Override
  public boolean safeFromWordTearing() {
    return true;
  }

  @Override
  public void save(ObjectOutputStream dos) throws IOException {
    dos.writeInt(IndexType.ARRAY.ordinal());
    dos.writeLong(mLength);
    dos.writeObject(mArray);
  }

  /**
   * Should only be called from {@link LongCreate#loadIndex(java.io.ObjectInputStream)}
   * @param ois stream to load from
   * @return index loaded from stream
   * @throws IOException if an IO error occurs
   */
  public static LongArray loadIndex(ObjectInputStream ois) throws IOException {
    final long length = ois.readLong();
    final long[] data;
    try {
      data = (long[]) ois.readObject();
    } catch (ClassNotFoundException e) {
      throw new IOException("Unrecognized index type: " + e.getMessage());
    }
    return new LongArray(data, length);
  }
}


