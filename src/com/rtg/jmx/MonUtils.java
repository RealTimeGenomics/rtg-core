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
package com.rtg.jmx;

import java.io.IOException;
import java.text.NumberFormat;

/**
 * Utility methods for formatting data from management beans.
 */
public final class MonUtils {

  private MonUtils() { }

  // Formatter with 0 fraction digits
  static final NumberFormat NF0 = NumberFormat.getInstance();
  // Formatter with 1 fraction digits
  static final NumberFormat NF1 = NumberFormat.getInstance();
  // Formatter with 2 fraction digits
  static final NumberFormat NF2 = NumberFormat.getInstance();
  static {
    NF0.setMinimumIntegerDigits(1);
    NF0.setMaximumFractionDigits(0);
    NF0.setMinimumFractionDigits(0);
    NF1.setMinimumIntegerDigits(1);
    NF1.setMaximumFractionDigits(1);
    NF1.setMinimumFractionDigits(1);
    NF2.setMinimumIntegerDigits(1);
    NF2.setMaximumFractionDigits(2);
    NF2.setMinimumFractionDigits(2);
  }

  static void padRight(Appendable out, String str, int width) throws IOException {
    padRight(out, str, width, ' ');
  }
  static void padRight(Appendable out, String str, int width, char pchar) throws IOException {
    out.append(str);
    for (int i = str.length(); i < width; i++) {
      out.append(pchar);
    }
  }

  static void pad(Appendable out, String str, int width) throws IOException {
    pad(out, str, width, ' ');
  }
  static void pad(Appendable out, String str, int width, char pchar) throws IOException {
    for (int i = str.length(); i < width; i++) {
        out.append(pchar);
    }
    out.append(str);
  }
}
