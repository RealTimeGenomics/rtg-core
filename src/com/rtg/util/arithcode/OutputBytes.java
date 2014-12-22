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

/**
 */
@TestClass("com.rtg.util.arithcode.BytesTest")
public final class OutputBytes implements Output {

  static final int BITS_PER_BYTE = 8;

  private final ByteChunks mOut;

  private int mByteBuffer;

  private int mBitsToGo;

  /**
   * @param out where the output bytes are written.
   */
  public OutputBytes(ByteChunks out) {
    mOut = out;
    reset();
  }

  private void flush() {
    mOut.append(mByteBuffer);
    reset();
  }

  private void reset() {
    mBitsToGo = BITS_PER_BYTE;
    mByteBuffer = 0;
  }

  private void writtenBit() {
    mBitsToGo--;
    if (mBitsToGo == 0) {
      flush();
    }
  }

  @Override
  public long endBlock() {
    if (mBitsToGo < BITS_PER_BYTE) {
      mByteBuffer = mByteBuffer << mBitsToGo;
      flush();
    }
    return mOut.length();
  }

  @Override
  public void close() {
    assert mBitsToGo == BITS_PER_BYTE;
    mOut.trim();
  }

  @Override
  public void writeBit(boolean bit) {
    if (bit) {
      writeBitTrue();
    } else {
      writeBitFalse();
    }
  }

  @Override
  public void writeBitTrue() {
    assert mBitsToGo > 0 && mBitsToGo <= BITS_PER_BYTE;
    mByteBuffer = (mByteBuffer << 1) | 1;
    writtenBit();
  }

  @Override
  public void writeBitFalse() {
    assert mBitsToGo > 0 && mBitsToGo <= BITS_PER_BYTE;
    mByteBuffer = mByteBuffer << 1;
    writtenBit();
  }
}
