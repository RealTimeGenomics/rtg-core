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

package com.rtg.variant.dna;

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.util.integrity.IntegerRange;

/**
 * Class for dealing with integer subranges representing
 * DNA values (<code>NACGTD</code>).
 * Uses an array to specify the single character values for each value.
 *
 */
@TestClass(value = {"com.rtg.variant.dna.DNARangeTest", "com.rtg.variant.dna.GeneralDNARangeTest"})
public class DNARange extends IntegerRange {

  private static final String VALUES = "NACGT";

  /** N unknown nucleotide. */
  public static final byte N = 0;

  /** A nucleotide. */
  public static final byte A = 1;

  /** C nucleotide. */
  public static final byte C = 2;

  /** G nucleotide. */
  public static final byte G = 3;

  /** T nucleotide. */
  public static final byte T = 4;

  /** Range for standard DNA encoding. */
  public static final DNARange RANGE = new DNARange(N, T);

  /**
   * Check that a nucleotide is in the correct range.
   * @param nt nucleotide to be checked.
   * @param lo first allowable value (inclusive).
   * @param hi last allowable value (inclusive).
   * @return true iff nt is in the allowed range.
   */
  public static boolean valid(final int nt, final int lo, final int hi) {
    assert N <= lo && lo <= hi && hi <= T;
    return lo <= nt && nt <= hi;
  }

  /**
   * @param nt nucleotide.
   * @return the complement of nt.
   */
  public static byte complement(final byte nt) {
    if (nt == N) {
      return N;
    }
    return (byte) (5 - nt);
  }


  /**
   * @param start first valid value (inclusive).
   * @param end last valid value (inclusive).
   */
  public DNARange(final int start, final int end) {
    super(false, -1, start, end);
  }


  @Override
  public String toString(int i) {
    return VALUES.substring(i, i + 1);
  }

  /**
   * Get a single character representation of nucleotide.
   * @param i value to be converted.
   * @return the character corresponding to the valid value i.
   */
  public char toChar(int i) {
    return VALUES.charAt(i);
  }

  @Override
  public int valueOf(String str) {
    final int res = VALUES.indexOf(str);
    if (res < 0) {
      throw new RuntimeException("Invalid string:" + str);
    }
    return res;
  }

  /**
   * Convert value of ch as an integer when interpreted as a DNA character.
   * Should obey the equation <code>i = valueOf(toString(i).charAt(0))</code>.
   * @param ch to be converted.
   * @return value.
   */
  public byte valueOf(char ch) {
    assert VALUES.length() <= Byte.MAX_VALUE;
    final byte res = (byte) VALUES.indexOf(ch);
    if (res < 0) {
      throw new RuntimeException("Invalid string:" + ch);
    }
    return res;
  }
}
