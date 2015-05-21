/*
 * Copyright (c) 2015. Real Time Genomics Limited.
 *
 * Use of this source code is bound by the Real Time Genomics Limited Software Licence Agreement
 * for Academic Non-commercial Research Purposes only.
 *
 * If you did not receive a license accompanying this file, a copy must first be obtained by email
 * from support@realtimegenomics.com.  On downloading, using and/or continuing to use this source
 * code you accept the terms of that license agreement and any amendments to those terms that may
 * be made from time to time by Real Time Genomics Limited.
 */
package com.rtg.reader;


/**
 * Utility functions for manipulating files that are not provided in the File
 * class.
 *
 */
public final class FastaUtils {

  /** ASCII PHRED highest allowed value */
  public static final int PHRED_UPPER_LIMIT_CHAR = '~';
  /** ASCII PHRED lowest allowed value */
  public static final int PHRED_LOWER_LIMIT_CHAR = '!';

  private FastaUtils() {
  }

  // Convert between raw and ascii phred quality values

  /**
   * @param quality the raw quality value
   * @return the ASCII phred value corresponding to raw quality value
   */
  public static char rawToAsciiQuality(byte quality) {
    return (char) (quality + PHRED_LOWER_LIMIT_CHAR);
  }

  /**
   * @param quality the ASCII quality value
   * @return the raw quality value corresponding to the ASCII quality value
   */
  public static byte asciiToRawQuality(char quality) {
    return (byte) (quality - PHRED_LOWER_LIMIT_CHAR);
  }

  /**
   * @param qualities the raw quality values
   * @return the ASCII phred values corresponding to raw quality values
   */
  public static char[] rawToAsciiQuality(byte[] qualities) {
    if (qualities == null) {
      return null;
    }
    final char[] result = new char[qualities.length];
    for (int i = 0; i < qualities.length; i++) {
      result[i] = rawToAsciiQuality(qualities[i]);
    }
    return result;
  }

  /**
   * Converts an array of bytes into a sanger-encoded quality string
   *
   * @param quality buffer containing input qualities
   * @return the quality string
   */
  public static String rawToAsciiString(byte[] quality) {
    return new String(rawToAsciiQuality(quality));
  }

  /**
   * Converts an array of bytes into a sanger-encoded quality string
   *
   * @param quality buffer containing input qualities
   * @param length the number of bytes from the input buffer to convert
   * @return the quality string
   */
  public static String rawToAsciiString(byte[] quality, int length) {
    final StringBuilder b = new StringBuilder();
    for (int i = 0; i < length; i++) {
      b.append(rawToAsciiQuality(quality[i]));
    }
    return b.toString();
  }

  /**
   * @param qualities the ASCII quality values
   * @return the raw quality values corresponding to the ASCII quality values
   */
  public static byte[] asciiToRawQuality(CharSequence qualities) {
    if (qualities == null) {
      return null;
    }
    final byte[] result = new byte[qualities.length()];
    for (int i = 0; i < result.length; i++) {
      result[i] = asciiToRawQuality(qualities.charAt(i));
    }
    return result;
  }

  /**
   * @param qualities the ASCII quality values
   * @return the raw quality values corresponding to the ASCII quality values
   */
  public static byte[] asciiToRawQuality(char[] qualities) {
    if (qualities == null) {
      return null;
    }
    final byte[] result = new byte[qualities.length];
    for (int i = 0; i < qualities.length; i++) {
      result[i] = asciiToRawQuality(qualities[i]);
    }
    return result;
  }

}
