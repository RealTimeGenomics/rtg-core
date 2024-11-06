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

package com.rtg.variant.bayes;

/**
 */
public class CodeHaploid implements Code {

  private final int mSize;

  /**
   * @param size number of haploid hypotheses.
   */
  public CodeHaploid(int size) {
    mSize = size;
  }

  @Override
  public int a(int n) {
    return n;
  }

  @Override
  public int b(int n) {
    throw new UnsupportedOperationException();
  }

  @Override
  public int bc(int n) {
    return n;
  }

  @Override
  public boolean homozygous(int n) {
    return true;
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
    return mSize;
  }

  @Override
  public int code(int a) {
    if (!valid(a)) {
      throw new IllegalArgumentException("a=" + a);
    }
    return a;
  }

  @Override
  public int code(int a, int b) {
    if (a == b) {
      return code(a);
    }
    throw new UnsupportedOperationException();
  }
}
