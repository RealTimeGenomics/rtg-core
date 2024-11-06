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

import java.io.IOException;

import com.rtg.util.integrity.IntegralAbstract;

/**
 */
public abstract class IntSet extends IntegralAbstract {

  private final IntSetCaller mCaller;

  /**
   * @param caller used to do calls for each unique value.
   */
  public IntSet(final IntSetCaller caller) {
    mCaller = caller;
  }

  /**
   * Add a new value to the set.
   * @param v the value to be added.
   * @throws IOException If an I/O error occurs
   */
  public abstract void add(final int v) throws IOException;

  /**
   * Call <code>call</code> for each unique value in the set.
   * Also remove all entries from the set.
   * @throws IOException If an I/O error occurs
   */
  public abstract void iterateClear() throws IOException;

  /**
   * Call <code>call</code> for each unique value in the set.
   * Also remove all entries from the set.
   * @throws IOException If an I/O error occurs
   */
  public abstract void iterateClearAll() throws IOException;

  /**
   * Called if the set overflows (exceeds capacity) or for each unique value
   * when <code>iterateClear</code> is called.
   * @param v value stored in set.
   * @throws IOException If an I/O error occurs
   */
  public void call(final int v) throws IOException {
    mCaller.call(v);
  }

}
