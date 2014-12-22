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
 * DNA values. Uses an array to specify the single character values for each
 * value.
 *
 */
@TestClass(value = {"com.rtg.variant.dna.DNARangeTest", "com.rtg.variant.dna.GeneralDNARangeTest"})
public class GeneralDNARange extends IntegerRange {

  private final String mValues;

  /**
   * @param values string of the (single) character values.
   * @param invalid true iff an invalid value is allowed.
   */
  GeneralDNARange(final String values, final boolean invalid) {
    super(invalid, -1, 0, values.length() - 1);
    mValues = values;
  }

  @Override
  public String toString(int i) {
    return mValues.substring(i, i + 1);
  }

  /**
   * Get a single character representation of nucleotide.
   * @param i value to be converted.
   * @return the character corresponding to the valid value i.
   */
  public char toChar(int i) {
    return mValues.charAt(i);
  }

  @Override
  public int valueOf(String str) {
    final int res = mValues.indexOf(str);
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
  public int valueOf(char ch) {
    final int res = mValues.indexOf(ch);
    if (res < 0) {
      throw new RuntimeException("Invalid string:" + ch);
    }
    return res;
  }
}
