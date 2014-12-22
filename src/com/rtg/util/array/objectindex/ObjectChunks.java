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

import com.rtg.util.integrity.Exam;

/**
 * Break array into chunks to fit within java convention that indices must
 * be ints.
 *
 * @param <A> type of objects stored in array.
 */
public class ObjectChunks<A> extends ObjectIndex<A> {

  private final int mChunkBits;

  private final int mChunkSize;

  private final int mChunkMask;

  private Object[][] mArray;

  /**
   * Constructs an index by splitting into array chunks.
   *
   * @param length of the index being created.
   */
  public ObjectChunks(final long length) {
    this(length, CHUNK_BITS);
  }

  /**
   * Constructs an index by splitting into array chunks.
   * This version sets the size of the chunks - it should only be used for testing.
   * @param length of the index being created.
   * @param chunkBits number of bits used for an entry in a chunk.
   */
  ObjectChunks(final long length, final int chunkBits) {
    super(length);
    assert chunkBits > 0 && chunkBits <= 30;
    mChunkBits = chunkBits;
    mChunkSize = 1 << mChunkBits;
    mChunkMask = mChunkSize - 1;

    final long ch = (length + mChunkSize - 1) / mChunkSize;
    if (ch > Integer.MAX_VALUE) {
      throw new RuntimeException("length requested too long length=" + length + " mChunkSize=" + mChunkSize);
    }
    final int chunks = (int) ch;
    mArray = new Object[chunks][];
    long left = mLength;
    for (int i = 0; i < chunks; i++) {
      final int assignedLength = left <= mChunkSize ? (int) left :  mChunkSize;
      if (assignedLength == 0) {
        throw new RuntimeException("zero assigned length");
      }
      mArray[i] = new Object[assignedLength];
      left -= assignedLength;
    }
    if (left != 0) {
      throw new RuntimeException("Did not assign requested memory mLength=" + mLength + " mChunkSize=" + mChunkSize + " left=" + left + " chunks=" + chunks);
    }
    assert integrity();
  }

  @Override
  public A get(final long index) {
    final int chunk = (int) (index >> mChunkBits);
    final int offset = (int) (index & mChunkMask);
    @SuppressWarnings("unchecked")
    final A ret = (A) mArray[chunk][offset];
    return ret;
  }

  @Override
  public void swap(final long index1, final long index2) {
    final int chunk1 = (int) (index1 >> mChunkBits);
    final int offset1 = (int) (index1 & mChunkMask);

    final int chunk2 = (int) (index2 >> mChunkBits);
    final int offset2 = (int) (index2 & mChunkMask);
    @SuppressWarnings("unchecked")
    final A temp = (A) mArray[chunk1][offset1];
    mArray[chunk1][offset1] = mArray[chunk2][offset2];
    mArray[chunk2][offset2] = temp;
  }

  @Override
  public void set(final long index, final A value) {
    final int chunk = (int) (index >> mChunkBits);
    final int offset = (int) (index & mChunkMask);
    mArray[chunk][offset] = value;
  }

  @Override
  public void close() {
    mArray = null;
  }

  /**
   * Get the chunk size. Should only be used for testing.
   *
   * @return the chunk size.
   */
  int chunkSize() {
    return mChunkSize;
  }

  @Override
  public boolean integrity() {
    super.integrity();
    Exam.assertTrue(mChunkBits > 0 && mChunkBits <= 31);
    Exam.assertTrue((mChunkSize & mChunkMask) == 0);
    Exam.assertTrue(mChunkSize > 0);
    Exam.assertTrue(mChunkMask > 0);
    Exam.assertTrue(mChunkMask + 1 == mChunkSize);
    if (mArray == null) { //close call has been made.
      return true;
    }
    long l = 0;
    for (Object[] aMArray : mArray) {
      l += aMArray.length;
    }
    Exam.assertTrue(l + ":" + mLength, l == mLength);
    return true;
  }
}


