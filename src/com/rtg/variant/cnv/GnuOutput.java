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
package com.rtg.variant.cnv;


import java.io.IOException;
import java.io.OutputStream;

import com.rtg.util.StringUtils;
import com.rtg.util.Utils;

/**
 */
public class GnuOutput implements LSMOutput, AutoCloseable {

  private static final String LS = StringUtils.LS;
  private final OutputStream mOut;

  /**
   * @param out output stream
   */
  public GnuOutput(final OutputStream out) {
    if (out == null) {
      throw new NullPointerException();
    }
    mOut = out;
  }

  protected static void writeValues(OutputStream out, final int start, int end, double sum, double sum2) throws IOException {
    //gnuplot
    //System.err.write(("gnuplot id=" + id + " start=" + start + " end=" + end + " sum=" + sum + " sum2=" + sum2 + LS).getBytes());
    final int length = end - start;
    final double mean = sum / length;
    final double threshold = sum * sum / length;
    final double stdDev = length == 1 ? 0.0 : Math.sqrt((sum2 - threshold) / (length - 1));
    //System.err.write(("gnuplot length=" + length + " mean=" + mean + LS).getBytes());

    final double stc = start - 1.25; //slightly left of (zero based) points
    final double endc = end - 1.75; //slightly right of points
    final String sts = Utils.realFormat(stc, 2);
    final String ends = Utils.realFormat(endc, 2);


    final double lo = mean - stdDev;
    final double hi = mean + stdDev;
    out.write((sts + " " + Utils.realFormat(mean, 3) + LS).getBytes());
    out.write((ends + " " + Utils.realFormat(mean, 3) + LS).getBytes());
    out.write(LS.getBytes());
    out.write((sts + " " + Utils.realFormat(lo, 3) + LS).getBytes());
    out.write((sts + " " + Utils.realFormat(hi, 3) + LS).getBytes());
    out.write(LS.getBytes());

    out.write((ends + " " + Utils.realFormat(lo, 3) + LS).getBytes());
    out.write((ends + " " + Utils.realFormat(hi, 3) + LS).getBytes());
    out.write(LS.getBytes());
  }

  @Override
  public void out(final String id, final int start, final int end, final double sum, final double sum2) throws IOException {
    writeValues(mOut, start, end, sum, sum2);
  }

  @Override
  public void close() throws IOException {
    mOut.close();
  }

  @Override
  public String toString() {
    return "GnuOutput";
  }
}
