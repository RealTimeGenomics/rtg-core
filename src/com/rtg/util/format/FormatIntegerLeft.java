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
package com.rtg.util.format;

/**
 * Format an integer, left justifying and padding it out to the
 * specified length as necessary. If it is too long leave it as is.
 * Optionally insert ","s. Usage: <code>FormatInteger format =
 * FormatInteger(3); .... System.out.println(format.format(-1);
 * //outputs " -1" .... StringBuilder sb = ... format.format(sb,-23);
 * //puts "-23" in string buffer.</code>
 *
 */
public class FormatIntegerLeft extends FormatInteger {

  /**
   * Construct a format for an integer.
   *
   * @param in number of positions (including - sign).
   */
  public FormatIntegerLeft(final int in) {
    this(in, false);
  }


  /**
   * Construct a format for a real.
   *
   * @param in number of positions before decimal point (including
   * - sign).
   * @param group true iff grouping to be used (as in 1,234 vs 1234)
   */
  public FormatIntegerLeft(final int in, final boolean group) {
    super(in, group);
  }


  /**
   * Format long into StringBuilder
   *
   * @param sb StringBuilder to which formatted value is to be appended
   * @param w long to be formatted
   * @return StringBuilder that has been modified (standard convention)
   */
  @Override
  public StringBuilder format(final StringBuilder sb, final long w) {
    final int initPosition = sb.length();
    sb.append(mLocalFormat.format(w));
    final int currLength = sb.length() - initPosition;

    //pad at right as necessary.
    final int pad = mLength - currLength;
    if (pad > 0) {
      sb.append(mPadding, 0, pad);
    }
    return sb;
  }
}
