/*
 * Copyright (c) 2014. Real Time Genomics Limited.
 *
 * Use of this source code is bound by the Real Time Genomics Limited Software Licence Agreement
 * for Academic Non-commercial Research Purposes only.
 *
 * If you did not receive a license accompanying this file, a copy must first be obtained by email
 * from support@realtimegenomics.com.  On downloading, using and/or continuing to use this source
 * code you accept the terms of that license agreement and any amendments to those terms that may
 * be made from time to time by Real Time Genomics Limited.
 */

package com.rtg.variant.sv;

import java.io.PrintStream;

import com.rtg.util.Utils;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;


/**
 */
public class SignalSum extends IntegralAbstract implements Signal {

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
    assert !Double.isNaN(sum) && !Double.isInfinite(sum);
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
  public boolean globalIntegrity() {
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
