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
 * Calculates the minimum or maximum complexity of a protein sequence.
 * Based on SEG algorithm. <code>K = 1 / L * log20 ( L! / prod(ni!))</code> for a sliding window of length L.
 *
 */
public class ProteinSequenceComplexity extends AbstractSequenceComplexity {

  // log_20(x!) for 0 to 20
  private static final double[] LOG_20_FACT = {0, 0, 0.231378, 0.598104, 1.06086, 1.5981, 2.19621, 2.84577, 3.5399, 4.27335, 5.04198, 5.84241, 6.6719, 7.5281, 8.40904, 9.31301, 10.2385, 11.1843, 12.1491, 13.132, 14.132};

  /**
   * Creates a protein sequence complexity calculator.
   *
   * @param length complexity region length
   */
  public ProteinSequenceComplexity(final int length) {
    super(length);
  }

  /**
   * Returns the complexity for the region of the sequence starting at <code>offset</code>.
   *
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

    final int[] baseCounts = new int[27]; // only needs to be 20 for proteins and 4 for nucleotides

    final int end = Math.min(sequence.length, offset + regionLength());
    final int length = end - offset;

    for (int i = offset; i < end; i++) {
      //System.err.println(i + " : " + sequence[i]);
      baseCounts[sequence[i]]++;
    }
    double logOfNFactorials = 0.0;
    for (int baseCount : baseCounts) {
      logOfNFactorials += LOG_20_FACT[baseCount];
    }
    //System.err.println("(" + LOG_20_FACT[length] + " - " + logOfNFactorials + ") / " + length);
    return (LOG_20_FACT[length] - logOfNFactorials) / length;
  }

  /**
   * Transform a human-readable Protein sequence into internal 0..26 bytes.
   * @return the encoded array
   */
  @Override
  protected byte[] encodeString(final String str) {
    final byte[] bytes = new byte[str.length()];
    for (int i = 0; i < str.length(); i++) {
      bytes[i] = (byte) str.charAt(i);
    }
    encodeArray(bytes);
    return bytes;
  }

  /**
   * Transform (in-situ) a human-readable Protein sequence into internal 0..26 bytes.
   * @return the encoded array (which will be the same array as <code>a</code>)
   */
  private byte[] encodeArray(final byte[] a) {
    for (int k = 0; k < a.length; k++) {
      final byte x = a[k];
      if (x >= (byte) 'a' && x <= (byte) 'z') {
        a[k] = (byte) (x - 'a' + 1);
      } else if (x >= (byte) 'A' && x <= (byte) 'Z') {
        a[k] = (byte) (x - 'A' + 1);
      } else {
        a[k] = 0;
      }
    }
    return a;
  }
}
