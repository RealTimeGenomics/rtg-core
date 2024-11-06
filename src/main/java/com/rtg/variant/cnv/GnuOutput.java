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
package com.rtg.variant.cnv;


import java.io.IOException;
import java.io.OutputStream;

import com.rtg.util.StringUtils;
import com.rtg.util.Utils;

/**
 */
public class GnuOutput implements LsmOutput, AutoCloseable {

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
