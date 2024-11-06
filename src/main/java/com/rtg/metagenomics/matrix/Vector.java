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
package com.rtg.metagenomics.matrix;

import java.util.Arrays;

import com.rtg.util.Utils;

/**
 */
public class Vector {

  private final int mSize;

  private final double[] mVector;

  /**
   * For testing.
   * @param v array of values to be converted to a vector.
   */
  public Vector(final double[] v) {
    mSize = v.length;
    mVector = v.clone();
  }

  /**
   * Copy constructor.
   * @param vec the Vector to copy
   */
  public Vector(final Vector vec) {
    mSize = vec.size();
    mVector = Arrays.copyOf(vec.mVector, mSize);
  }

  /**
   * @param size length of vector.
   */
  public Vector(final int size) {
    mSize = size;
    mVector = new double[size];
  }

  /**
   * Get the <code>i'th</code> value.
   * @param i  index.
   * @return the value.
   */
  public double get(final int i) {
    return mVector[i];
  }

  /**
   * Set the <code>i'th</code> value.
   * @param i index.
   * @param v value.
   */
  public void set(final int i, final double v) {
    mVector[i] = v;
  }

  /**
   * Increment the <code>i'th</code> value.
   * @param i index.
   * @param v value.
   */
  public void incr(final int i, final double v) {
    mVector[i] += v;
  }

  /**
   * Get the length of the vector.
   * @return the length (&ge; 0).
   */
  public int size() {
    return mSize;
  }

  @Override
  public String toString() {
    final StringBuilder sb = new StringBuilder();
    for (int i = 0; i < mSize; ++i) {
      sb.append("  ").append(Utils.realFormat(mVector[i], 4));
    }
    return sb.toString();
  }
}
