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

import java.io.PrintStream;

import com.rtg.util.Utils;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;


/**
 */
public final class SignalSum extends IntegralAbstract implements Signal {

  private final Signal[] mSignals;

  private final PrintStream mOut;

  private final String mColumnName;


  /**
   * @param out a <code>PrintStream</code> to which distribution contributions will be written.
   * @param column the name of the column for output.
   * @param records array of <code>SAMRecord</code> objects used to construct the signal.
   */
  public SignalSum(PrintStream out, String column, Signal... records) {
    mSignals = records;
    mOut = out;
    mColumnName = column;
    assert globalIntegrity();
  }

  /**
   * @param column the name of the column for output.
   * @param records array of <code>SAMRecord</code> objects used to construct the signal.
   */
  public SignalSum(String column, Signal... records) {
    mOut = null;
    mSignals = records;
    mColumnName = column;
    assert globalIntegrity();
  }

  @Override
  public double value(int position) {
    double sum = 0;
    if (mOut != null) {
      mOut.print(position + "\t");
    }
    for (Signal mSignal : mSignals) {
      final double contrib = mSignal.value(position);
      sum += contrib;
      if (mOut != null) {
        mOut.print(Utils.realFormat(contrib, 2) + "\t");
      }
    }
    if (mOut != null) {
      mOut.println(Utils.realFormat(sum, 2));
    }
    assert Double.isFinite(sum);
    return sum;
  }

  @Override
  public void toString(StringBuilder sb) {
    sb.append("Sum:");
  }

  @Override
  public boolean integrity() {
    Exam.assertNotNull(mSignals);
    return true;
  }

  @Override
  public final boolean globalIntegrity() {
    integrity();
    for (Signal mSignal : mSignals) {
      Exam.globalIntegrity(mSignal);
    }
    return true;
  }

  @Override
  public String columnLabel() {
    return mColumnName;
  }

}
