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

package com.rtg.util.arithcode;


import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.util.array.byteindex.ByteChunks;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;

/**
 */
@TestClass("com.rtg.util.arithcode.BytesTest")
public final class InputBytes extends IntegralAbstract implements Input {

  private final ByteChunks mIn;

  private final long mStart;

  private final long mEnd;

  private long mPosition;

  private int mBitsToGo = 0;

  private int mByteBuffer;

  private boolean mEof = false;

  /**
   * @param in index containing concatenated compressed blocks.
   * @param start index of first byte in current block (0 based).
   * @param end index one past last byte in current block (0 based).
   */
  public InputBytes(ByteChunks in, long start, long end) {
    mIn = in;
    mStart = start;
    mEnd = end;
    mPosition = start;
    fetchInput();
  }

  @Override
  public boolean readBit() {
    if (mEof) {
      return false;
    }
    mBitsToGo--;
    final boolean res = (mByteBuffer & (1 << mBitsToGo)) != 0;
    fetchInput();
    return res;
  }

  private void fetchInput() {
    if (mBitsToGo == 0) {
      if (mPosition == mEnd) {
        mEof = true;
      } else {
        mByteBuffer = mIn.getByte(mPosition);
        mBitsToGo = OutputBytes.BITS_PER_BYTE;
        mPosition++;
      }
    }
  }

  @Override
  public boolean integrity() {
    Exam.assertTrue(toString(), 0 <= mStart && mStart <= mPosition && mPosition <= mEnd && mEnd <= mIn.length());
    Exam.assertTrue(0 <= mBitsToGo && mBitsToGo <= OutputBytes.BITS_PER_BYTE);
    Exam.assertTrue(0 <= mByteBuffer && mByteBuffer < (1 << OutputBytes.BITS_PER_BYTE));
    return true;
  }

}
