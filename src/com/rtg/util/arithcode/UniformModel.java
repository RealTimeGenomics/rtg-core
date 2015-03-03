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

/**
 * A singleton uniform distribution byte model.  Provides a single
 * static member that is a non-adaptive model assigning equal
 * likelihood to all 256 bytes.
 */
public final class UniformModel implements DetailedModel {

  /**
   * A re-usable uniform model.
   */
  public static final UniformModel MODEL = new UniformModel(256);

  /**
   * Construct a uniform model.
   * @param numOutcomes the number of alternatives to encode
   */
  UniformModel(final int numOutcomes) {
    mNumOutcomes = numOutcomes;
  }

  private final int mNumOutcomes;

  @Override
  public int totalCount() {
    return mNumOutcomes;
  }

  @Override
  public void encode(ArithEncoder encoder, int symbol) {
    encoder.encode(symbol, symbol + 1, mNumOutcomes);
  }

  @Override
  public int decode(ArithDecoder decoder) {
    final int symbol = decoder.getCurrentSymbolCount(mNumOutcomes);
    decoder.removeSymbolFromStream(symbol, symbol + 1, mNumOutcomes);
    return symbol;
  }

  @Override
  public int pointToSymbol(int midCount) {
    return midCount;
  }

  @Override
  public void interval(int symbol, int[] result) {
    result[0] = symbol;
    result[1] = symbol + 1;
    result[2] = mNumOutcomes;
  }
}
