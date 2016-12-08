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
package com.rtg.complexity;

/**
 * Abstract class for SEG complexity calculation classes.
 *
 */
public abstract class AbstractSequenceComplexity implements SequenceComplexity {
  // maximum region for complexity calculations
  protected static final int MAX_LENGTH = 20;

  // region length for complexity calculations
  private final int mRegionLength;

  /**
   * Constructs a complexity calculator for the specified region length.
   *
   * @param regionLength length of region for complexity calculation
   */
  public AbstractSequenceComplexity(int regionLength) {
    if (regionLength < 1 || regionLength > MAX_LENGTH) {
      throw new IllegalArgumentException("Complex region length must be between 1 and 20: " + regionLength);
    }
    mRegionLength = regionLength;
  }

  /**
   * Returns the complexity for the region of the sequence starting at <code>offset</code>.
   *
   * @param sequence sequence to operate on
   * @param offset point to calculate from
   * @return complexity for the sequence
   */
  protected abstract double complexity(byte[] sequence, int offset);

  /**
   * Returns a byte encoded version of the given sequence.
   *
   * @param sequence sequence to encode
   * @return encoded sequence
   */
  protected abstract byte[] encodeString(final String sequence);


  /**
   * Return complexity region length.
   *
   * @return region length
   */
  protected int regionLength() {
    return mRegionLength;
  }


  /**
   * Returns the minimum complexity of any region in the <code>sequence</code>.
   *
   * @param sequence sequence to operate on
   * @return minimum complexity
   */
  @Override
  public final double minComplexity(byte[] sequence) {
    double res = -1.0;
    for (int i = 0; i < Math.max(1, sequence.length - mRegionLength); ++i) {
      final double c = complexity(sequence, i);
      if (res == -1.0 || c < res) {
        res = c;
      }
    }
    return res;
  }

  /**
   * Returns the minimum complexity of any region in the <code>sequence</code>.
   *
   * @param sequence sequence to operate on
   * @return minimum complexity
   */
  @Override
  public double maxComplexity(String sequence) {
    return maxComplexity(encodeString(sequence));
  }

  /**
   * Returns the maximum complexity of any region in the <code>sequence</code>.
   *
   * @param sequence sequence to operate on
   * @return minimum complexity
   */
  @Override
  public final double maxComplexity(byte[] sequence) {
    double res = -1.0;
    for (int i = 0; i < Math.max(1, sequence.length - mRegionLength); ++i) {
      final double c = complexity(sequence, i);
      if (c > res) {
        res = c;
      }
    }
    return res;
  }

  /**
   * Returns the maximum complexity of any region in the <code>sequence</code>.
   *
   * @param sequence sequence to operate on
   * @return minimum complexity
   */
  @Override
  public double minComplexity(String sequence) {
    return minComplexity(encodeString(sequence));
  }

}
