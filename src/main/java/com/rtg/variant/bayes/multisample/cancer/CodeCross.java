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

package com.rtg.variant.bayes.multisample.cancer;

import com.rtg.variant.bayes.Code;

/**
 * A non-commutative code which is a Cartesian cross product of another code.
 */
public class CodeCross implements Code {

  private final int mBaseSize;
  private final int mSize;

  /**
   * @param size of the underlying code whose product is being taken.
   */
  CodeCross(final int size) {
    assert size > 0 && size <= 46340;
    mBaseSize = size;
    mSize = mBaseSize * mBaseSize;
  }

  @Override
  public int a(int n) {
    return n / mBaseSize;
  }

  @Override
  public int b(int n) {
    return n % mBaseSize;
  }

  @Override
  public int bc(int n) {
    return b(n);
  }

  @Override
  public boolean homozygous(int n) {
    return (n / mBaseSize) == (n % mBaseSize);
  }

  @Override
  public boolean valid(int hyp) {
    return 0 <= hyp && hyp < mSize;
  }

  @Override
  public int size() {
    return mSize;
  }

  @Override
  public int rangeSize() {
    return mBaseSize;
  }

  @Override
  public int code(int a) {
    if (a >= mBaseSize) {
      throw new IllegalArgumentException(String.valueOf(a));
    }
    return a;
  }

  @Override
  public int code(int a, int b) {
    return a * mBaseSize + b;
  }
}
