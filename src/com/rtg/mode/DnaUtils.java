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
package com.rtg.mode;



/**
 * DNA utility functions.
 *
 */
public final class DnaUtils {

  private static final char[] BASES_LOWER = {
    Character.toLowerCase(DNA.values()[0].toString().charAt(0)),
    Character.toLowerCase(DNA.values()[1].toString().charAt(0)),
    Character.toLowerCase(DNA.values()[2].toString().charAt(0)),
    Character.toLowerCase(DNA.values()[3].toString().charAt(0)),
    Character.toLowerCase(DNA.values()[4].toString().charAt(0)),
  };
  private static final char[] BASES = {
    Character.toUpperCase(DNA.values()[0].toString().charAt(0)),
    Character.toUpperCase(DNA.values()[1].toString().charAt(0)),
    Character.toUpperCase(DNA.values()[2].toString().charAt(0)),
    Character.toUpperCase(DNA.values()[3].toString().charAt(0)),
    Character.toUpperCase(DNA.values()[4].toString().charAt(0)),
  };

  /** Value for the unknown residue */
  public static final byte UNKNOWN_RESIDUE = 0;

  private DnaUtils() { }

  /**
   * Convert a number into a string of nucleotide codes.
   * This is handy for displaying hash function bits.
   * @param x number to be displayed
   * @param len number of bits and codes to display.
   * @return string giving bit decomposition of x.
   */
  public static String toCodes(final long x, final int len) {
    if (len <= 0 || len > 32) {
      throw new IllegalArgumentException("length out of range=" + len);
    }
    final StringBuilder sb = new StringBuilder();
    final int left = 64 - 2 * len;
    long t = x << left;
    long u = x << (left + len);
    for (int i = 0; i < len; i++) {
      final int code = (t < 0 ? 2 : 0) + (u < 0 ? 1 : 0);
      sb.append(new String[] {"A", "C", "G", "T"}[code]);
      t = t << 1;
      u = u << 1;
    }
    assert u == 0;
    return sb.toString();
  }

  /**
   * Convert a number into a string of nucleotide codes.
   * This is handy for displaying hash function bits.
   * @param v0 first set of bits.
   * @param v1 second set of bits.
   * @param len number of bits and codes to display.
   * @return string giving bit decomposition of x.
   */
  public static String toCodes(final long v0, final long v1, final int len) {
    if (len <= 0 || len > 64) {
      throw new IllegalArgumentException("length out of range=" + len);
    }
    final StringBuilder sb = new StringBuilder();
    final int left = 64 - len;
    long t = v0 << left;
    long u = v1 << left;
    for (int i = 0; i < len; i++) {
      final int code = (t < 0 ? 2 : 0) + (u < 0 ? 1 : 0);
      sb.append(new String[] {"A", "C", "G", "T"}[code]);
      t = t << 1;
      u = u << 1;
    }
    assert u == 0;
    return sb.toString();
  }

  /**
   * Take reverse complement of a string.
   * @param seq to be reverse complemented.
   * @return the reversed sequence.
   */
  public static String reverseComplement(final String seq) {
    final StringBuilder sb = new StringBuilder();
    for (int i = seq.length() - 1; i >= 0; i--) {
      final char c = seq.charAt(i);
      switch (c) {
      case 'a':
        sb.append('t');
        break;
      case 'A':
        sb.append('T');
        break;
      case 'c':
        sb.append('g');
        break;
      case 'C':
        sb.append('G');
        break;
      case 'g':
        sb.append('c');
        break;
      case 'G':
        sb.append('C');
        break;
      case 't':
        sb.append('a');
        break;
      case 'T':
        sb.append('A');
        break;
      case '-':
        break;
      default: // n
        sb.append(c);
        break;
      }
    }
    return sb.toString();
  }

  /**
   * Convert from a 0-4 based numerical base to an uppercase <code>NACGT</code> style character.
   * @param b byte representing the numerical value of the base
   * @return character representation of the byte (<code>NACGT</code>)
   */
  public static char getBase(final int b) {
    return BASES[b];
  }

  /**
   * Converts an encoded base into a uppercase <code>ACGTN</code> character.
   * Handles all the off-the-end cases by returning 'N'.
   *
   * @param a array of encoded residues.
   * @param p position of the desired base in <code>a</code>.
   * @return a character form of the base, or 'N'.
   */
  public static char base(final byte[] a, final int p) {
    return p >= 0 && p < a.length ? getBase(a[p]) : getBase(UNKNOWN_RESIDUE);
  }

  /**
   * Convert from a 0-4 based numerical base to a lowercase <code>nacgt</code> style character.
   * Note, we don't generally output lowercase bases, so if possible use <code>getBase</code> instead.
   * @param b byte representing the numerical value of the base
   * @return character representation of the byte (<code>nacgt</code>)
   */
  public static char getBaseLower(final int b) {
    return BASES_LOWER[b];
  }

  /**
   * Convert a binary DNA sequence to a human readable string using uppercase characters
   * @param seq DNA in internal 0..4 bytes
   * @param start offset into array
   * @param length length of DNA
   * @return readable string
   */
  public static String bytesToSequenceIncCG(final byte[] seq, final int start, final int length) {
    final StringBuilder sb = new StringBuilder();
    for (int i = start; i < start + length; i++) {
      if (i < 0 || i >= seq.length) {
        sb.append('N');
      } else if (seq[i] != 5) { //allow for CG spacer, but strip it
        sb.append(BASES[seq[i]]);
      }
    }
    return sb.toString();
  }

  /**
   * Convert a binary DNA sequence to a human readable string using uppercase characters
   * @param seq DNA in internal 0..4 bytes
   * @return readable string
   */
  public static String bytesToSequenceIncCG(final byte[] seq) {
    return bytesToSequenceIncCG(seq, 0, seq.length);
  }

  /**
   * Transform a human-readable DNA sequence into internal 0..4 bytes.
   * @param str Eg. <code>"ACGTN"</code> will become {1,2,3,4,0}.
   * @return the encoded array
   */
  public static byte[] encodeString(final String str) {
    return DnaUtils.encodeArray(str.getBytes());
  }

  /**
   * Transform a human-readable DNA sequence into internal 0..4 bytes.
   * Hyphens are ignored and can be used to space sequence data for readability in tests.
   * @param str Eg. <code>"ACGTN"</code> will become {1,2,3,4,0}.
   * @return the encoded array
   */
  public static byte[] encodeStringWithHyphen(final String str) {
    return DnaUtils.encodeString(str.replace("-", ""));
  }

  /**
   * Transform (in-situ) a human-readable DNA sequence into internal 0..4 bytes.
   * @param a Eg. {'a','c','g','t','n'} will become {1,2,3,4,0}.
   * @param length length to convert
   * @return the encoded array (which will be the same array as <code>a</code>)
   */
  public static byte[] encodeArray(final byte[] a, final int length) {
    for (int k = 0; k < length; k++) {
      switch (a[k]) {
      case (byte) 'a':
      case (byte) 'A':
        a[k] = 1;
        break;
      case (byte) 'c':
      case (byte) 'C':
        a[k] = 2;
        break;
      case (byte) 'g':
      case (byte) 'G':
        a[k] = 3;
        break;
      case (byte) 't':
      case (byte) 'T':
        a[k] = 4;
        break;
      default:
        a[k] = 0;
        break;
      }
    }
    return a;
  }

  /**
   * Transform a human-readable DNA sequence into the reverse complemented version, uppercase.
   * @param src Eg. {'a','c','g','t','n'} will become {'N', 'A', 'C', 'G', 'T'}.
   * @param dest destination byte array
   * @param length length to convert
   */
  public static void reverseComplement(byte[] src, byte[] dest, int length) {
    for (int k = 0; k < length; k++) {
      switch (src[k]) {
      case (byte) 'a':
      case (byte) 'A':
        dest[length - 1 - k] = 'T';
        break;
      case (byte) 'c':
      case (byte) 'C':
        dest[length - 1 - k] = 'G';
        break;
      case (byte) 'g':
      case (byte) 'G':
        dest[length - 1 - k] = 'C';
        break;
      case (byte) 't':
      case (byte) 'T':
        dest[length - 1 - k] = 'A';
        break;
      default:
        dest[length - 1 - k] = 'N';
        break;
      }
    }
  }

  /**
   * Transform (in-situ) a human-readable DNA sequence into internal 0..4 bytes.
   * @param a Eg. {'a','c','g','t','n'} will become {1,2,3,4,0}.
   * @return the encoded array (which will be the same array as <code>a</code>)
   */
  public static byte[] encodeArray(final byte[] a) {
    return encodeArray(a, a.length);
  }

  /**
   * Convert a FASTQ string to a byte array of phred scores.
   * @param string the quality string to convert
   * @return byte array of phred scores
   */
  public static byte[] fastqToPhred(final String string) {
    final byte[] bytes = new byte[string.length()];
    for (int i = 0; i < string.length(); i++) {
      bytes[i] = (byte) (string.charAt(i) - 33);
    }
    return bytes;
  }

  /**
   * Convert an array of bytes, in which each byte is a binary phred quality
   * score, to printable ASCII representation of the quality scores.
   * @param buffer Array of bytes in which each byte is a binary phred score.
   * @param offset Where in buffer to start conversion.
   * @param length How many bytes of buffer to convert.
   * @param reverse if the binary phred score is to be reversed
   * @return String with ASCII representation of those quality scores.
   */
  public static String phredToFastq(final byte[] buffer, final int offset, final int length, final boolean reverse) {
    final char[] chars = new char[length];
    for (int i = 0; i < length; i++) {
      chars[reverse ? length - i - 1 : i] = (char) ((buffer[offset + i]) + 33);
    }
    return new String(chars);
  }

}
