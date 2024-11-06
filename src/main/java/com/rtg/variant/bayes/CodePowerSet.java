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
 * Power set without empty set.
 */
public class CodePowerSet implements Code {

  private final int mSize;

  /**
   * Code representing the power set of specified number of items but excluding empty set.
   * @param n underlying description size
   */
  public CodePowerSet(final int n) {
    assert n <= 31; // Actually n > 6 or so will make for slow calling ...
    mSize = (1 << n) - 1;
  }

  @Override
  public int a(final int n) {
    throw new UnsupportedOperationException();
  }

  @Override
  public int b(final int n) {
    throw new UnsupportedOperationException();
  }

  @Override
  public int bc(final int n) {
    throw new UnsupportedOperationException();
  }

  @Override
  public boolean homozygous(final int n) {
    return ((n + 1) & n) == 0;
  }

  @Override
  public boolean valid(final int hyp) {
    return 0 <= hyp && hyp < mSize;
  }

  @Override
  public int size() {
    return mSize;
  }

  @Override
  public int rangeSize() {
    throw new UnsupportedOperationException();
  }

  @Override
  public int code(final int a) {
    throw new UnsupportedOperationException();
  }

  @Override
  public int code(final int a, final int b) {
    throw new UnsupportedOperationException();
  }
}
