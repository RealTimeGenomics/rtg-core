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
package com.rtg.util.array.byteindex;

/**
 * Index is implemented using a straight forward long array.
 * This is so that short instances of ShortIndex can be as efficient as possible.
 *
 */
public final class ByteArray extends ByteIndex {

  private final byte[] mArray;

  /**
   * This should be called directly only in tests.
   *
   * @param length number of items to be stored.
   */
  public ByteArray(final long length) {
    super(length);
    //assert length <= MAX_LENGTH;
    mArray = new byte[(int) length];
  }

  @Override
  public byte getByte(final long index) {
    final int ii = (int) index;
    if (ii != index) {
      throw new IndexOutOfBoundsException(String.valueOf(index));
    }
    return mArray[ii];
  }

  @Override
  public void setByte(final long index, final byte value) {
    final int ii = (int) index;
    if (ii != index) {
      throw new IndexOutOfBoundsException(String.valueOf(index));
    }
    mArray[ii] = value;
  }

  @Override
  public boolean integrity() {
    super.integrity();
    assert mArray.length == mLength;
    return true;
  }

  @Override
  public boolean safeFromWordTearing() {
    return true;
  }
}


