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

import java.util.Arrays;

/**
 * <P>Model where counts for each symbol are provided externally.
 *
 * @version 1.0
 */
public final class StaticModel implements DetailedModel {

  private final int mTotalCount;

  private final int[] mCounts;

  @Override
  public int totalCount() {
    return mTotalCount;
  }

  @Override
  public int pointToSymbol(int midCount) {
    final int bs = Arrays.binarySearch(mCounts, midCount);
    if (bs >= 0) {
      return bs;
    }
    return -bs - 2;
  }

  @Override
  public void interval(int symbol, int[] result) {
    result[0] = mCounts[symbol];
    result[1] = mCounts[symbol + 1];
    result[2] = mTotalCount;
  }

  @Override
  public void encode(ArithEncoder encoder, int symbol) {
    encoder.encode(mCounts[symbol], mCounts[symbol + 1], mTotalCount);
  }

  @Override
  public int decode(ArithDecoder decoder) {
    final int cnt = decoder.getCurrentSymbolCount(mTotalCount);
    final int symbol = pointToSymbol(cnt);
    decoder.removeSymbolFromStream(mCounts[symbol], mCounts[symbol + 1], mTotalCount);
    return symbol;
  }

  /**
   * Construct a uniform model.
   */
  StaticModel(final int[] counts) {
    mCounts = new int[counts.length + 1];
    int sofar = 0;
    for (int i = 1; i <= counts.length; i++) {
      sofar += counts[i - 1] + 1;
      mCounts[i] = sofar;
    }
    mTotalCount = mCounts[mCounts.length - 1];
  }

}
