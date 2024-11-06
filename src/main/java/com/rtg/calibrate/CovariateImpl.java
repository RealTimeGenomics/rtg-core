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

package com.rtg.calibrate;

import com.reeltwo.jumble.annotations.TestClass;

/**
 */
@TestClass({"com.rtg.calibrate.CovariateReadGroupTest", "com.rtg.calibrate.CovariateBaseQualityTest"})
public abstract class CovariateImpl implements Covariate {

  private final String mName;

  private int mSize;
  protected int mNewSize;

  /**
   * Default constructor for <code>CovariateImpl</code>
   * @param name the name of the variable (lower case)
   * @param size the initial size of the variable
   */
  public CovariateImpl(String name, int size) {
    assert size >= 0;
    mSize = mNewSize = size;
    mName = name;
  }

  @Override
  public String name() {
    return mName;
  }

  @Override
  public final int size() {
    return mSize;
  }

  @Override
  public void resized() {
    mSize = newSize();
  }

  @Override
  public final int newSize() {
    return mNewSize;
  }

  protected final void setNewSize(int size) {
    assert size >= mNewSize;
    mNewSize = size;
  }

  @Override
  public final boolean sizeChanged() {
    assert newSize() >= mSize;
    return newSize() > mSize;
  }

  /**
   * Default implementation returns number as string.
   */
  @Override
  public String valueString(int value) {
    return Integer.toString(value);
  }

  /**
   * Default implementation returns string parsed as number.
   */
  @Override
  public int parse(String value) {
    final int ret = Integer.parseInt(value);
    if (ret >= newSize()) {
      setNewSize(ret + 1);
    }
    return ret;
  }

  /**
   * Default implementation returns string parsed as number.
   */
  @Override
  public int valueOf(String value) {
    return parse(value);
  }

  @Override
  public String toString() {
    return name() + ":" + size() + ":" + newSize();
  }


}
