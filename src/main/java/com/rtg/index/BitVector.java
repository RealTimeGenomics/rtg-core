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
/**
 * A vector of bits.
 */
public final class BitVector extends AbstractBitVector {

  private final int[] mArray;

  /**
   * Creates a new <code>BitVector</code> with the specified length.
   *
   * @param length the number of bits.
   */
  public BitVector(final long length) {
    super(length);
    final long allocLength = (length + MASK) >> BITS_PER_ENTRY;
    if (allocLength >= (Integer.MAX_VALUE >> 2)) {
      throw new RuntimeException("BitVector is too long " + length);
    }
    mArray = new int[(int) allocLength];
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
    final int v = mArray[(int) x];
    final int s = 1 << (int) (index & MASK);
    return (v & s) != 0;
  }

  @Override
  public void set(final long index) {
    checkBounds(index);
    final long x = index >> BITS_PER_ENTRY;
    mArray[(int) x] |= 1 << (int) (index & MASK);
  }

  @Override
  public void reset(final long index) {
    checkBounds(index);
    final long x = index >> BITS_PER_ENTRY;
    mArray[(int) x] &= ~(1 << (int) (index & MASK));
  }

  @Override
  public long bytes() {
    return (((long) mArray.length) << BITS_PER_ENTRY) >> BITS_PER_BYTE;
  }
}


