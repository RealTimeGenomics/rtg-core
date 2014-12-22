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
package com.rtg.index.hash.ngs.general;


import static com.rtg.util.StringUtils.TAB;

import java.io.IOException;
import java.io.PrintStream;

import com.rtg.util.format.FormatReal;

/**
 */
public final class MaskAnalysis {

  private static final double SUBSTITUTION_PROB = 0.02;

  private static final FormatReal FMT = new FormatReal(7, 2);

  private static final FormatReal PERC = new FormatReal(2, 2);

  private MaskAnalysis() { }

  /**
   * @param args command line arguments: read length, substitutions, indels, indel length
   * @throws IOException If an I/O error occurs
   */
  public static void main(final String[] args) throws IOException {
    final int r = Integer.parseInt(args[0]);
    final int s = Integer.parseInt(args[1]);
    final int i = Integer.parseInt(args[2]);
    final int l = Integer.parseInt(args[3]);
    final PrintStream out = System.out;
    out.println("read=" + r + " substitutions=" + s + " indels=" + i + " indelLength=" + l);
    out.println("   w" + TAB + "build" + TAB + "    search" + TAB + "search1" + TAB + "search2" + TAB + "    miss");
    for (int w = 1; w <= Math.min(32, r - s); w++) {
      final Skeleton sk = new Skeleton(r, w, s, i, l);
      if (!sk.valid()) {
        continue;
      }
      //System.err.println(sk);
      final long b = sk.numberWindows();
      final int c1 = MaskIndelCount.indelCount(sk);
      final int rl = sk.windowLength();
      final double c2 = Math.pow(3.5, -rl) * 1.0e9;
      final double q = Math.exp(-sk.chunkLength() * SUBSTITUTION_PROB);
      final double p = 1.0 - q;
      double miss = 0.0;
      final int n = r / sk.chunkLength();
      for (int j = s + 1; j <= n; j++) {
        miss += Math.pow(q, n - j) * Math.pow(p, j) * Util.binomial(n, j);
      }
      out.println(w + TAB + b + TAB + FMT.format(c1 + c2) + TAB + c1 + TAB + FMT.format(c2) + TAB + PERC.format(miss * 100.0) + "%");
    }
  }
}
