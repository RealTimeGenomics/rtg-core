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
package com.rtg.index;

import com.rtg.util.array.intindex.IntChunks;


/**
 * A vector of bits.
 * Uses underlying big arrays so that can be longer than 2^34 bits.
 */
public final class BigBitVector extends AbstractBitVector {

  private final IntChunks mArray;

  /**
   * Creates a new <code>BitVector</code> with the specified length.
   *
   * @param length the number of bits.
   */
  public BigBitVector(final long length) {
    super(length);
    final long allocLength = (length + MASK) >> BITS_PER_ENTRY;
    mArray = new IntChunks(allocLength);
  }

  private void checkBounds(final long index) {
    if (index >= mLength || index < 0) {
      throw new ArrayIndexOutOfBoundsException(index + ":" + mLength);
    }
  }

  @Override
  public boolean get(final long index) {
    checkBounds(index);
    final long x = index >> BITS_PER_ENTRY;
    final int v = mArray.getInt(x);
    final int s = 1 << (int) (index & MASK);
    return (v & s) != 0;
  }

  @Override
  public void set(final long index) {
    checkBounds(index);
    final long x = index >> BITS_PER_ENTRY;
    final int v = mArray.getInt(x) | 1 << (int) (index & MASK);
    mArray.setInt(x, v);
  }

  @Override
  public void reset(final long index) {
    checkBounds(index);
    final long x = index >> BITS_PER_ENTRY;
    final int v = mArray.getInt(x) & ~(1 << (int) (index & MASK));
    mArray.setInt(x, v);
  }


  @Override
  public long bytes() {
    return (mArray.length() << BITS_PER_ENTRY) >> BITS_PER_BYTE;
  }
}


