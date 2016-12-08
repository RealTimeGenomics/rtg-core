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

import com.rtg.mode.DnaUtils;

/**
 * Calculates the minimum or maximum complexity of a DNA sequence.
 * Based on SEG algorithm. <code>K = 1 / L * log4 ( L! / prod(ni!))</code> for a sliding window of length L.
 *
 */
public class DnaSequenceComplexity extends AbstractSequenceComplexity {
  // log_4(x!) for 0 to 20
  private static final double[] LOG_4_FACT = {0, 0, 0.5, 1.29248, 2.29248, 3.45345, 4.74593, 6.1496, 7.6496, 9.23457, 10.8955, 12.6252, 14.4177, 16.2679, 18.1716, 20.1251, 22.1251, 24.1688, 26.2538, 28.3777, 30.5387};

  /**
   * Creates a DNA sequence complexity calculator.
   *
   * @param length complexity region length
   */
  public DnaSequenceComplexity(final int length) {
    super(length);
  }

  /**
   * Returns the complexity for the region of the sequence starting at <code>offset</code>.
   * The sequence is assumed to have been encoded as returned by <code>Utils.encodeArray()</code>
   *
   * @param sequence sequence to operate on
   * @param offset point to calculate from
   * @return complexity for the sequence
   */
  @Override
  protected double complexity(byte[] sequence, int offset) {
    if (sequence == null || sequence.length == 0 || offset >= sequence.length) {
      return -1.0;
    }

    final int[] baseCounts = new int[6]; // a,c,g,t,n,space

    final int end = Math.min(sequence.length, offset + regionLength());
    final int length = end - offset;

    for (int i = offset; i < end; ++i) {
      baseCounts[sequence[i]]++;
    }
    double logOfNFactorials = 0.0;
    for (int baseCount : baseCounts) {
      //System.err.println(i + " : " + baseCounts[i]);
      logOfNFactorials += LOG_4_FACT[baseCount];
    }
    //System.err.println("(" + LOG_4_FACT[length] + " - " + logOfNFactorials + ") / " + length);
    return (LOG_4_FACT[length] - logOfNFactorials) / length;
  }

  @Override
  protected byte[] encodeString(String sequence) {
    return DnaUtils.encodeString(sequence);
  }

}
