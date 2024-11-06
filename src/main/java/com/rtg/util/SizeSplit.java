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
package com.rtg.util;

/**
 * Facilities for splitting a range as evenly as possible.
 * Handy for doing things like allocating work to threads.
 */
public class SizeSplit {
  private final int mN;

  private final int mD;

  private final int mS;

  private final int mR;

  /**
   * @param n number of unit to be divided.
   * @param d the number of groups they are to be divided into.
   */
  public SizeSplit(final int n, final int d) {
    assert n >= 0 && d >= 1;
    mN = n;
    mD = d;
    mS = n / d;
    mR = n - d * mS;
    assert mR >= 0;
    assert mS >= 0;
  }

  /**
   * Get the first position (zero based, inclusive) allocated to "unit" i.
   * @param i unit of work ranging from 0 to d inclusive.
   * @return the last allocated value (zero based, inclusive).
   */
  public int start(final int i) {
    assert 0 <= i && i <= mD;
    return i <= mR ? i * (mS + 1) : mR + i * mS;
  }

  @Override
  public String toString() {
    return "SizeSplit: n=" + mN + " d=" + mD + " s=" + mS + " r=" + mR;
  }

}
