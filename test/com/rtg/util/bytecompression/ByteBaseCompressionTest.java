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

import com.rtg.util.array.longindex.LongChunks;

/**
 */
public class ByteBaseCompressionTest extends AbstractByteCompressionTest {

  @Override
  protected ByteCompression getCompressor() {
    return new ByteBaseCompression(10);
  }

  public void testNoCompression() {
    check(new ByteBaseCompression(256), "1234", "094123", "86754", "", "", "12", "1");
  }

  public void testBytes() {
    final ByteBaseCompression cmp = new ByteBaseCompression(256);
    final String[] strings = {"1234", "094123", "86754", "", "", "12", "1"};
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
      final int length = pointers[i + 1] - pointers[i];
      cmp.add(bytes, pointers[i], length);
      assertEquals(length, cmp.length(i));
    }
    assertEquals(82, cmp.bytes());
    cmp.freeze();
    assertEquals(82, cmp.bytes());
    for (int i = 0; i < pointers.length - 1; i++) {
      final int length = pointers[i + 1] - pointers[i];
      assertEquals(length, cmp.length(i));
    }
  }

  public void testSdfConstructor() {
    final String[] strings = {"1234", "094123", "86754", "", "", "12", "1"};
    final LongChunks pointers = new LongChunks(strings.length + 1);
    pointers.set(0, 0);
    int numBytes = 0;
    for (int i = 0; i < strings.length; i++) {
      final int length = strings[i].length();
      numBytes += length;
      pointers.set(i + 1, numBytes);
    }
    final ByteArray data = new BitwiseByteArray(numBytes, 6);
    for (int i = 0; i < strings.length; i++) {
      final long start = pointers.get(i);
      for (int j = 0; j < strings[i].length(); j++) {
        data.set(start + j, new byte[] {(byte) (strings[i].charAt(j) - '0')}, 1);
      }
    }
    final ByteCompression cmp = new ByteBaseCompression(data, pointers);
    try {
      cmp.add(new byte[0], 0, 0);
      fail();
    } catch (final RuntimeException e) {
      assertEquals("Adding to a frozen ByteCompression", e.getMessage());
    }
    final long totalBytes = cmp.bytes();
    cmp.freeze();
    assertTrue(totalBytes >= cmp.bytes());
    for (int i = 0; i < pointers.length() - 1; i++) {
      final long start = pointers.get(i);
      final int length = (int) (pointers.get(i + 1) - start);
      final byte[] outBytes = new byte[length];
      cmp.get(outBytes, i, 0, length);
      for (int j = 0; j < length; j++) {
        assertEquals(data.get(j + start), outBytes[j]);
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
