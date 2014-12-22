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
 * <P>Performs arithmetic decoding, converting bit input into
 * cumulative probability interval output.  Returns probabilities as
 * integer counts <code>low</code>, <code>high</code> and
 * <code>total</code>, with the range being
 * <code>[low/total,high/total)</code>.
 *
 * <P>For more details, see <a href="../../../tutorial.html">The Arithmetic Coding Tutorial</a>.
 *
 * @version 1.1
 * @see ArithEncoder
 * @see Input
 */
@TestClass("com.rtg.util.arithcode.ArithTest")
public final class ArithDecoder extends ArithCoder {

  /**
   * Construct an arithmetic decoder that reads from the given
   * bit input.
   * @param in Bit input from which to read bits.
   */
  public ArithDecoder(Input in) {
    mIn = in;
    for (int i = 1; i <= CODE_VALUE_BITS; ++i) {
      bufferBit();
    }
  }

  /**
   * Returns a count for the current symbol that will be between the
   * low and high counts for the symbol in the model given the total count.
   * Once symbol is retrieved, the model is used to compute the actual low,
   * high and total counts and {@link #removeSymbolFromStream} is called.
   * @param totalCount The current total count for the model.
   * @return A count that is in the range above or equal to the low count and less than the high count of the next symbol decoded.
   */
  public int getCurrentSymbolCount(int totalCount) {
    return (int) (((mValue - mLow + 1) * totalCount - 1) / (mHigh - mLow + 1));
  }

  /**
   * Removes a symbol from the input stream.  Called after {@link #getCurrentSymbolCount}.
   * @param lowCount Cumulative count for symbols indexed below symbol to be removed.
   * @param highCount <code>lowCount</code> plus count for this symbol.
   * @param totalCount Total count for all symbols seen.
   */
  public void removeSymbolFromStream(long lowCount, long highCount, long totalCount) {
    final long range = mHigh - mLow + 1;
    mHigh = mLow + (range * highCount) / totalCount - 1;
    mLow = mLow + (range * lowCount) / totalCount;
    while (true) {
      if (mHigh < HALF) {
        // no effect
      } else if (mLow >= HALF) {
        mValue -= HALF;
        mLow -= HALF;
        mHigh -= HALF;
      } else if (mLow >= FIRST_QUARTER && mHigh <= THIRD_QUARTER) {
        mValue -= FIRST_QUARTER;
        mLow -= FIRST_QUARTER;
        mHigh -= FIRST_QUARTER;
      } else {
        return;
      }
      mLow <<= 1; // = 2 * mLow;
      mHigh = (mHigh << 1) + 1; // 2 * mHigh + 1;
      bufferBit();
    }
  }

  /**
   * Input stream from which to read bits.
   */
  private final Input mIn;

  /**
   * Current bits for decoding.
   */
  private long mValue; // implied = 0;

  /**
   * Reads a bit from the underlying bit input stream and buffers it.
   */
  private void bufferBit() {
    mValue = mValue << 1;
    if (mIn.readBit()) {
      ++mValue;
    }
  }

}
