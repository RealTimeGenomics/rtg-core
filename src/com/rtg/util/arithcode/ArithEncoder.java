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


/**
 * <P>Performs arithmetic encoding, converting cumulative probability
 * interval input into bit output.  Cumulative probability intervals
 * are given as integer counts <code>low</code>, <code>high</code> and
 * <code>total</code>, with the range being
 * <code>[low/total,high/total)</code>.
 *
 * <P>For more details, see <a href="../../../tutorial.html">The Arithmetic Coding Tutorial</a>.
 * </P>
 * <P>
 * The public methods below can be called in the following order:
 * <br>
 * <code>(endBlock* encode* endBlock+)* close endBlock*</code>
 * </P>
 * @version 1.1
 * @see ArithDecoder
 * @see Output
 */
@TestClass("com.rtg.util.arithcode.ArithTest")
public final class ArithEncoder extends ArithCoder {

  /**
   * Construct an arithmetic coder from a bit output.
   * @param out Underlying bit output.
   */
  public ArithEncoder(Output out) {
    mOut = out;
  }

  /**
   * Get a key which is passed into Input to retrieve a particular block.
   * @return the key associated with the next block.
   */
  public long endBlock() {
    if (mState == State.ENCODING) {
      ++mBitsToFollow; // need a final bit (not sure why)
      if (mLow < FIRST_QUARTER) {
        bitPlusFollowFalse();
      } else {
        bitPlusFollowTrue();
      }
      mBitsToFollow = 0;
      mLow = 0;
      mHigh = TOP_VALUE;
      mState = State.END_BLOCK;
    }
    return mOut.endBlock();
  }

  /**
   * Close the arithmetic encoder, writing all bits that are
   * buffered and closing the underlying output streams.
   */
  public void close() {
    assert mState == State.END_BLOCK;
    mState = State.CLOSED;
    mOut.close();
  }

  /**
   * Encodes an interval expressed as a low count, high count and total count.
   * The high count is taken to be exclusive, and the resulting range is
   * <code>highCount - lowCount + 1</code>.
   * @param lowCount Cumulative count of symbols below current one.
   * @param highCount Cumulative count of symbols below current one plus current one.
   * @param totalCount Cumulative count of all symbols.
   */
  public void encode(int lowCount, int highCount, int totalCount) {
    assert mState == State.ENCODING || mState == State.END_BLOCK;
    mState = State.ENCODING;
    final long range = mHigh - mLow + 1;
    mHigh = mLow + (range * highCount) / totalCount - 1;
    mLow  = mLow + (range * lowCount) / totalCount;
    while (true) {
      if (mHigh < HALF) {
        bitPlusFollowFalse();
      } else if (mLow >= HALF) {
        bitPlusFollowTrue();
        mLow -= HALF;
        mHigh -= HALF;
      } else if (mLow >= FIRST_QUARTER && mHigh < THIRD_QUARTER) {
        ++mBitsToFollow;
        mLow -= FIRST_QUARTER;
        mHigh -= FIRST_QUARTER;
      } else {
        return;
      }
      mLow <<= 1;
      mHigh = (mHigh << 1) + 1;
    }
  }

  /**
   * Bit output stream for writing encoding bits.
   */
  private final Output mOut;

  /**
   * Number of bits beyond first bit that were normalized.
   */
  private int mBitsToFollow; // implied = 0;

  private static enum State {
    CLOSED,
    END_BLOCK,
    ENCODING
  }

  private State mState = State.END_BLOCK;

  /**
   * Write a <code>true</code> bit, and then a number of <code>false</code> bits
   * equal to the number of bits to follow.
   */
  private void bitPlusFollowTrue() {
    for (mOut.writeBitTrue(); mBitsToFollow > 0; --mBitsToFollow) {
      mOut.writeBitFalse();
    }
  }

  /**
   * Write a <code>false</code> bit, and then a number of <code>true</code> bits
   * equal to the number of bits to follow.
   */
  private void bitPlusFollowFalse() {
    for (mOut.writeBitFalse(); mBitsToFollow > 0; --mBitsToFollow) {
      mOut.writeBitTrue();
    }
  }

}
