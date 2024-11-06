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
package com.rtg.variant.sv;

import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;

/**
 */
public class SignalCount extends IntegralAbstract implements Signal {

  private final SamCounts mRecords;

  private final int mWindowLo;

  private final int mWindowHi;

  private final String mColumnName;

  /**
   * @param records array of <code>SAMRecord</code> objects used to construct the signal.
   * @param lo low index of the convolution window.
   * @param hi  high index of the convolution window.
   * @param column the name of the column for output.
   */
  public SignalCount(SamCounts records, int lo, int hi, String column) {
    mRecords = records;
    mWindowLo = lo;
    mWindowHi = hi;
    mColumnName = column;
    assert integrity();
  }

  @Override
  public double value(int position) {
    return mRecords.count(position, mWindowLo, mWindowHi);
  }

  @Override
  public void toString(StringBuilder sb) {
    sb.append("Count:").append(mWindowLo).append(":").append(mWindowHi);
  }

  @Override
  public final boolean integrity() {
    Exam.assertNotNull(mRecords);
    Exam.assertTrue(mWindowLo < 0);
    Exam.assertTrue(mWindowHi > 0);
    return true;
  }

  @Override
  public String columnLabel() {
    return mColumnName;
  }


}
