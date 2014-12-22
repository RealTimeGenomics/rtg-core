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
  public boolean integrity() {
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
