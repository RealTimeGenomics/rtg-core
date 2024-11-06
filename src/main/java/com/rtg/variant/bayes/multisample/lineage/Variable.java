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
package com.rtg.variant.bayes.multisample.lineage;

import com.rtg.util.Utils;

/**
 * Name of a discrete random variable.
 */
public class Variable implements Comparable<Variable> {

  private final int mSize;
  private final String mName;

  /**
   * Construct a new discrete random variable with a name and size.
   * @param name name of variable
   * @param size number of hypotheses for this variable
   */
  Variable(final String name, final int size) {
    mSize = size;
    mName = name;
  }

  /**
   * Number of entries in the factor.
   * @return size of factor
   */
  int size() {
    return mSize;
  }

  @Override
  public boolean equals(Object that) {
    if (!(that instanceof Variable)) {
      return false;
    }
    final Variable s = (Variable) that;
    return mSize == s.mSize && mName.equals(s.mName);
  }

  @Override
  public int hashCode() {
    return Utils.pairHash(mName.hashCode(), mSize);
  }

  @Override
  public String toString() {
    return mName;
  }

  @Override
  public int compareTo(Variable o) {
    return mName.compareTo(o.mName);
  }
}
