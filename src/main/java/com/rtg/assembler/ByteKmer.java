/*
 * Copyright (c) 2018. Real Time Genomics Limited.
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

package com.rtg.assembler;

import com.rtg.assembler.graph.Contig;
import com.rtg.mode.DNA;

/**
 * Represent a kmer as an array of bit packed bytes.
 * Each byte can contain 4 bases.
 * The final 2 bits of the last byte in the array indicate how many bases are stored in the last byte (0-3)
 */
public class ByteKmer extends AbstractKmer {
  static final int CHARS_PER_BYTE = Byte.SIZE / 2;
  byte[] mKmer;


  ByteKmer(byte[] b) {
    mKmer = b;
  }

  ByteKmer(String kmer) {
    mKmer = toPacked(kmer);
  }

  /**
   * Pack bases from a string into our byte[] representation
   * @param s string representation of the bases
   * @return a byte array packed according to the <code>ByteKmer</code> scheme
   */
  static byte[] toPacked(String s) {
    final byte[] packed = array(s.length());
    for (int i = 0; i < s.length(); ++i) {
      final int base = DNA.valueOf(s.charAt(i)).ordinal() - 1 ;
      setPos(packed, i, base);
    }
    setSizeBits(packed, s.length());
    return packed;
  }

  /**
   * Pack bases from a 1 byte per base array into our compact byte[] representation
   * @param b 1 byte per base representation of the bases
   * @param start first position in the array to scan
   * @param end end position, exclusive
   * @return a byte array packed according to the <code>ByteKmer</code> scheme
   */
  static byte[] toPacked(byte[] b, int start, int end) {
    //    System.err.println(Arrays.toString(b) + " " +  start + "-" + end);
    final int size = end - start;
    final byte[] packed = array(size);
    //    System.err.println(packed.length);
    for (int i = start; i < end; ++i) {
      final int base = b[i] - 1 ;
      setPos(packed, i - start, base);
    }
    setSizeBits(packed, size);
    return packed;
  }
  static byte[] toPacked(Contig contig, int start, int end) {
    final int size = end - start;
    final byte[] packed = array(size);
    for (int i = start; i < end; ++i) {
      final int base = contig.nt(i) - 1;
      setPos(packed, i - start, base);
    }
    setSizeBits(packed, size);
    return packed;
  }

  /**
   * Create an array large enough to hold a kmer of the specified size
   * @param size the kmer size
   */
  private static byte[] array(double size) {
    return new byte[(int) Math.ceil((size + 1) / CHARS_PER_BYTE)];
  }

  /**
   * @return the size of the kmer in bases
   */
  @Override
  public int length() {
    final byte lastByte = mKmer[mKmer.length - 1];
    return (mKmer.length - 1) * CHARS_PER_BYTE + ((lastByte >> 6) & 3);
  }

  /**
   * set the base at position <code>pos</code> to <code>val</code> in the packed byte[] <code>b</code>
   * @param b the byte[] to alter
   * @param pos the base index to change
   * @param val the new value
   */
  static void setPos(byte[] b, int pos, int val) {
    if (pos < 0 || b.length * CHARS_PER_BYTE <= pos) {
      throw new RuntimeException("pos=" + pos);
    }
    final int byteIndex = pos / CHARS_PER_BYTE;
    final int bitPos = pos % CHARS_PER_BYTE;
    byte mask = (byte) 0xff;
    mask = (byte) (mask ^ (((byte) 3) << (bitPos * 2)));
    b[byteIndex] = (byte) (b[byteIndex] & mask);
    b[byteIndex] = (byte) ((byte) (val << (bitPos * 2)) | b[byteIndex]);
  }

  /**
   * set the bit which indicate the size of the kmer
   * @param b the byte[] to alter
   * @param size length of the kmer
   */
  static void setSizeBits(byte[] b, int size) {
    final byte lastByte = b[b.length - 1];
    // Move the bit size into the right place
    final byte sizeBits = (byte) ((size & 3) << 6);
    // zero out the previous size bits
    final byte ntBits = (byte) (((1 << 6) - 1) & lastByte);
    b[b.length - 1] = (byte) (ntBits + sizeBits);

  }

  /**
   * Retrieve the base at position <code>pos</code> in the packed byte[] <code>b</code>
   * @param b a byte array packed according to our scheme
   * @param pos base index to retrieve
   * @return the base as an int
   */
  static int getPos(byte[] b, int pos) {
    if (pos < 0 || b.length * CHARS_PER_BYTE <= pos) {
      throw new RuntimeException("pos=" + pos);
    }
    final int byteIndex = pos / CHARS_PER_BYTE;
    final int bitPos = pos % CHARS_PER_BYTE;
    final byte val = (byte) (b[byteIndex] >> (bitPos * 2));
    return val & (byte) 3;
  }

  @Override
  public Kmer successor(byte nt) {
    final byte base = (byte) (nt - 1);
    final int size = length();
    final byte[] newBytes = new byte[mKmer.length];
    setPos(newBytes, size - 1, base);
    for (int i = 0; i < size - 1; ++i) {
      setPos(newBytes, i, baseAt(i + 1));
    }
    setSizeBits(newBytes, size);
    return new ByteKmer(newBytes);
  }

  @Override
  public Kmer predecessor(byte nt) {
    final byte base = (byte) (nt - 1);
    final int size = length();
    final byte[] newBytes = new byte[mKmer.length];
    for (int i = 1; i < size; ++i) {
      setPos(newBytes, i, baseAt(i - 1));
    }
    setPos(newBytes, 0, base);
    setSizeBits(newBytes, size);
    return new ByteKmer(newBytes);
  }

  static int complement(int base) {
    return 3 - base;
  }

  @Override
  public Kmer minimalKmer() {
    boolean isMinimal = true;
    for (int i = 0; i < length(); ++i) {
      final int val = baseAt(i);
      final int reverse = complement(baseAt(length() - i - 1));
      if (val < reverse) {
        break;
      }
      if (val > reverse) {
        isMinimal = false;
        break;
      }
    }
    if (isMinimal) {
      return this;
    } else {
      return reverse();
    }
  }

  @Override
  public Kmer reverse() {
    final byte[] b = new byte[mKmer.length];
    final int size = length();
    for (int i = 0; i < size; ++i) {
      setPos(b, i, 3 - getPos(mKmer, size - i - 1));
    }
    setSizeBits(b, size);
    return new ByteKmer(b);
  }

  int baseAt(int pos) {
    return getPos(mKmer, pos);
  }

  @Override
  public byte nt(int pos) {
    return (byte) (getPos(mKmer, pos) + 1);
  }

  private static class Factory implements KmerFactory {
    @Override
    public Kmer make(byte[] kmer, int start, int end) {
      return new ByteKmer(toPacked(kmer, start, end));
    }
    @Override
    public Kmer make(Contig contig, int start, int end) {
      return new ByteKmer(toPacked(contig, start, end));
    }

  }

  /**
   * @return a factory that will construct <code>ByteKmer</code>s from strings
   */
  static KmerFactory factory() {
    //Diagnostic.userLog("Using byte kmer factory");
    return new Factory();
  }
}
