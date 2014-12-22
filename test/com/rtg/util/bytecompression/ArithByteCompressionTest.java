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

import com.rtg.util.arithcode.Order0ModelBuilder;

/**
 */
public class ArithByteCompressionTest extends AbstractByteCompressionTest {

  @Override
  protected ByteCompression getCompressor() {
    return new ArithByteCompression(10, 20, new Order0ModelBuilder(10));
  }

  public void testBytes() {
    final ArithByteCompression cmp = new ArithByteCompression(10, 12, new Order0ModelBuilder(10));
    cmp.integrity();
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
      cmp.integrity();
    }
    assertEquals(74, cmp.bytes());
    cmp.freeze();
    cmp.integrity();
    assertEquals(74, cmp.bytes());
  }

  //set initial count high enough that counts never frozen
  public void testBytesUncompressed() {
    final ArithByteCompression cmp = new ArithByteCompression(10, 1000, new Order0ModelBuilder(10));
    cmp.integrity();
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
      cmp.integrity();
    }
    assertEquals(80, cmp.bytes());
    cmp.freeze();
    cmp.integrity();
    assertEquals(74, cmp.bytes());
  }


}
