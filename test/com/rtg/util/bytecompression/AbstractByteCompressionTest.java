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

package com.rtg.util.bytecompression;

import java.util.Random;

import junit.framework.TestCase;

/**
 */
public abstract class AbstractByteCompressionTest extends TestCase {

  /**
   * Any compressor returned by this function should build its
   * compression tables by the time 10 byte arrays.
   */
  protected abstract ByteCompression getCompressor();

  protected static String[] randomStrings(final long seed, final int number, final int maxLength) {
    final Random randomLen = new Random(seed);
    final Random randomCh = new Random(13 * seed + 101);

    final String[] strs = new String[number];
    for (int i = 0; i < number; i++) {
      final int length = randomLen.nextInt(maxLength);
      final StringBuilder sb = new StringBuilder();
      for (int j = 0; j < length; j++) {
        final int ch = randomCh.nextInt(10);
        sb.append(ch);
      }
      strs[i] = sb.toString();
    }
    return strs;
  }

  public void test() {
    check(getCompressor(), "", "012345678", "1234412", "12");
    check(getCompressor(), "012345678", "1234412", "12", "");
    check(getCompressor(), "012345678", "1234412", "12", "", "123", "352", "31", "563", "1337", "7", "11");
    check(getCompressor(), randomStrings(42, 20, 100));
  }

  /**
   * Check that the same strings that go in come out again.
   * The input strings to this method will have '0' subtracted from them to convert
   * to bytes.
   * @param cmp the compressor to use
   * @param strings the strings to compress
   */
  protected static void check(ByteCompression cmp, String... strings) {
    final int[] pointers = new int[strings.length + 1];
    pointers[0] = 0;
    int numBytes = 0;
    for (int i = 0; i < strings.length; i++) {
      final int length = strings[i].length();
      numBytes += length;
      pointers[i + 1] = numBytes;
    }
    final byte[] bytes = new byte[numBytes];
    for (int i = 0; i < strings.length; i++) {
      final int start = pointers[i];
      for (int j = 0; j < strings[i].length(); j++) {
        bytes[start + j] = (byte) (strings[i].charAt(j) - '0');
      }
    }
    for (int i = 0; i < pointers.length - 1; i++) {
      cmp.add(bytes, pointers[i], pointers[i + 1] - pointers[i]);
    }
    final long totalBytes = cmp.bytes();
    cmp.freeze();
    assertTrue(totalBytes >= cmp.bytes());
    for (int i = 0; i < pointers.length - 1; i++) {
      final int start = pointers[i];
      final int length = pointers[i + 1] - start;
      final byte[] outBytes = new byte[length];
      cmp.get(outBytes, i, 0, length);
      for (int j = 0; j < length; j++) {
        assertEquals(bytes[j + start], outBytes[j]);
      }
    }
    //System.err.println(cmp.toString());
    for (int i = 0; i < pointers.length - 1; i++) {
      final int start = pointers[i];
      final int length = pointers[i + 1] - start;
      for (int j = 0; j < length; j++) {
        //System.err.println("i=" + i + " j=" + j);
        final byte[] outByte = new byte[1];
        cmp.get(outByte, i, j, 1);
        assertEquals(bytes[j + start], outByte[0]);
      }
    }
    try {
      cmp.add(new byte[0], 0, 0);
      fail();
    } catch (final RuntimeException e) {
      assertEquals("Adding to a frozen ByteCompression", e.getMessage());
    }
  }
}
