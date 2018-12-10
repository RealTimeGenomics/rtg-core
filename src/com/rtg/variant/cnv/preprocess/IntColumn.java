/*
 * Copyright (c) 2016. Real Time Genomics Limited.
 *
 * Use of this source code is bound by the Real Time Genomics Limited Software Licence Agreement
 * for Academic Non-commercial Research Purposes only.
 *
 * If you did not receive a license accompanying this file, a copy must first be obtained by email
 * from support@realtimegenomics.com.  On downloading, using and/or continuing to use this source
 * code you accept the terms of that license agreement and any amendments to those terms that may
 * be made from time to time by Real Time Genomics Limited.
 */
package com.rtg.variant.cnv.preprocess;

/**
 * Convenience wrapper for data to be parsed / formatted as integer
 */
class IntColumn extends NumericColumn {

  IntColumn(String name) {
    super(name, "%d");
  }

  @Override
  String toString(int i) {
    return String.format(mFormat, (int) get(i));
  }

  @Override
  protected double toDouble(String strValue) {
    double value;
    try {
      value = Integer.parseInt(strValue);
    } catch (NumberFormatException e) {
      value = Double.NaN;
    }
    return value;
  }
}
