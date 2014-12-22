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
package com.rtg.util.array;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.List;

import com.rtg.util.array.longindex.LongIndex;
import com.rtg.util.io.ByteArrayIOUtils;
import com.rtg.util.io.FileUtils;

/**
 * Array utility functions.
 *
 */
public final class ArrayUtils {

  private ArrayUtils() { }

  /**
   * Sum an array of longs.
   *
   * @param a array
   * @return sum
   */
  public static long sum(final long[] a) {
    long sum = 0;
    for (final long x : a) {
      sum += x;
    }
    return sum;
  }

  /**
   * Sum an array of ints.  Returns a long in order to minimise changes of overflow.
   *
   * @param a array
   * @return sum
   */
  public static long sum(final int[] a) {
    long sum = 0;
    for (final int x : a) {
      sum += x;
    }
    return sum;
  }

  /**
   * Read an array of longs from a file.
   *
   * @param f file
   * @return array
   * @exception IOException if an I/O error occurs
   */
  public static long[] readLongArray(final File f) throws IOException {
    final byte[] b = new byte[(int) f.length()];
    try (InputStream stream = new FileInputStream(f)) {
      final int length = stream.read(b);
      if (length != b.length) {
        throw new IOException();
      }
    }
    return ByteArrayIOUtils.convertToLongArray(b);
  }

  /**
   * Read an array of longs from a file.
   *
   * @param f file
   * @param a index to read into
   * @param offset in the array.
   * @param addend amount to add to each entry
   * @return length
   * @exception IOException if an I/O error occurs
   */
  public static int readInts(final File f, final LongIndex a, final long offset, final long addend) throws IOException {
    return readInts(f, 0, (int) f.length() / 4, a, offset, addend);
  }

  /**
   * Read an array of longs from a file.
   *
   * @param f file
   * @param fileOffset offset into the file to start reading at
   * @param fileEnd position in the file to finish
   * @param a index to read into
   * @param offset in the array.
   * @param addend amount to add to each entry
   * @return length
   * @exception IOException if an I/O error occurs
   */
  public static int readInts(final File f, final int fileOffset, final int fileEnd, final CommonIndex a, final long offset, final long addend) throws IOException {
    final int length = (fileEnd - fileOffset) * 4;
    final byte[] b = new byte[length];
    try (InputStream stream = new FileInputStream(f)) {

      // Skip has failed lets attempt to use read to fix things up.
      int remaining = fileOffset * 4;
      long skipped;
      while (remaining > 0 && (skipped = FileUtils.streamSkip(stream, remaining)) > 0) {
        remaining -= (int) skipped;
      }
      if (remaining > 0) {
        throw new IOException();
      }
      int soFar = 0;
      int read;
      while (soFar < length && (read = stream.read(b, soFar, length - soFar)) > 0) {
        soFar += read;
      }
      if (soFar < length) {
        throw new IOException();
      }
    }
    for (int i = 0; i < b.length; i += 4) {
      a.set(offset + i / 4, ByteArrayIOUtils.bytesToIntBigEndian(b, i) + addend);
    }
    return b.length / 4;
  }

  /**
   * Return list of boxed bytes as a primitive array.
   *
   * @param l list
   * @return byte array
   */
  public static byte[] asByteArray(final List<Byte> l) {
    final byte[] a = new byte[l.size()];
    for (int i = 0; i < a.length; i++) {
      a[i] = l.get(i);
    }
    return a;
  }

  /**
   * Return list of boxed longs as a primitive array.
   *
   * @param l list
   * @return long array
   */
  public static long[] asLongArray(final List<Long> l) {
    final long[] a = new long[l.size()];
    for (int i = 0; i < a.length; i++) {
      a[i] = l.get(i);
    }
    return a;
  }

  /**
   * @param input list of integers represented as decimals separated by comma.
   * @return the values separated into an array
   * @throws NumberFormatException if the string format is invalid
   */
  public static int[] parseIntArray(final String input) {
    final String[] strings = input.split(", *");
    final int[] values = new int[strings.length];
    for (int i = 0; i < strings.length; i++) {
      values[i] = Integer.parseInt(strings[i]);
    }
    return values;
  }

  /**
   * Reverse a byte array in place.
   * @param arr the array to reverse.
   */
  public static void reverseArrayInPlace(final byte[] arr) {
    for (int i = 0, j = arr.length - 1; i < j; i++, j--) {
      final byte temp = arr[i];
      arr[i] = arr[j];
      arr[j] = temp;
    }
  }

}
