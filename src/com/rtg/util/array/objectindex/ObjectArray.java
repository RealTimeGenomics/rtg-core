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
package com.rtg.util.array.objectindex;

/**
 * Index is implemented using a straight forward long array.
 * This is so that short instances of IntIndex can be as efficient as possible.
 *
 */
class ObjectArray<A> extends ObjectIndex<A> {

  private Object[] mArray;

  ObjectArray(final long length) {
    super(length);
    assert length <= MAX_LENGTH;
    mArray = new Object[(int) length];
  }

  @Override
  public A get(final long index) {
    final int ii = (int) index;
    if (ii != index) {
      throw new IndexOutOfBoundsException(String.valueOf(index));
    }
    @SuppressWarnings("unchecked")
    final A ret = (A) mArray[ii];
    return ret;
  }

  @Override
  public void set(final long index, final A value) {
    final int ii = (int) index;
    if (ii != index) {
      throw new IndexOutOfBoundsException(String.valueOf(index));
    }
    mArray[ii] = value;
  }

  @Override
  public void close() {
    mArray = null;
  }

  @Override
  public boolean integrity() {
    super.integrity();
    assert mArray == null || mArray.length == mLength;
    return true;
  }
}


