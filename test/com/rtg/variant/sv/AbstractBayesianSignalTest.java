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

import static com.rtg.util.StringUtils.TAB;

import java.io.PrintStream;

import com.rtg.sam.ReadGroupUtils;
import com.rtg.util.Utils;
import com.rtg.util.diagnostic.Diagnostic;

import junit.framework.TestCase;

/**
 * Common code for testing <code>BayesianSignal</code>s.
 * Includes generating <code>SamArray</code>s with distributions that should be exactly findable.
 */
public abstract class AbstractBayesianSignalTest extends TestCase {

  protected static final double PROPER_RATE = 1.0;
  protected static final int READ_LENGTH = 20;
  protected static final double FRAGMENT_MEAN = 65.0;
  protected static final double FRAGMENT_STD_DEV = 10.0;
  protected static final double GAP_MEAN = 25.0;
  protected static final double GAP_STD_DEV = 10.0;
  protected static final int MAX_ALIGNMENT = 5;
  protected static final double PROPER_RANDOM_RATE = 0.01;
  protected static final double DISCORDANT_RATE = 0.01;
  protected static final int WIDTH = 3;
  protected static final double UNMATED_RATE = 0.001;
  protected static final ReadGroupStats STATS = new ReadGroupStats(ReadGroupUtils.UNKNOWN_RG, READ_LENGTH, FRAGMENT_MEAN, FRAGMENT_STD_DEV, GAP_MEAN, GAP_STD_DEV, MAX_ALIGNMENT, PROPER_RATE, PROPER_RANDOM_RATE, DISCORDANT_RATE, UNMATED_RATE);
  protected static final int WINDOW_LO = STATS.lo();
  protected static final int WINDOW_HI = STATS.hi();
  protected static final int WINDOW = WINDOW_HI - WINDOW_LO;
  protected static final int BREAK = WINDOW - WINDOW_LO;

  protected static int rev(int x, boolean reverse) {
    if (!reverse) {
      return x;
    }
    final int rev = -x - READ_LENGTH;
    assert WINDOW_LO <= rev && rev < WINDOW_HI : "x=" + x + " reverse=" + reverse + " rev=" + rev + " READ_LENGTH=" + READ_LENGTH + " WINDOW_LO=" + WINDOW_LO + " WINDOW_HI=" + WINDOW_HI;
    if (rev < WINDOW_LO) {
      return WINDOW_LO;
    }
    //System.err.println("x=" + x + " reverse=" + reverse + " offset=" + rev);
    return rev;
  }

  @Override
  public void setUp() throws Exception {
    Diagnostic.setLogStream();
  }

  SamArray constant(final double v) {
    final SamArray sa = new SamArray(2 * BREAK);
    for (int i = 0; i < 2 * BREAK; i++) {
      sa.increment(i, v);
    }
    return sa;
  }

  public void testBreakDetection() {
    check(false/*left BreakPoint*/, false);
    check(true/*right BreakPoint*/, false);
  }

  private static final int BREAK_TOLERANCE = 1;

  private void check(final boolean reverse, final boolean debug) {
    checkSignal(localSignal(reverse, debug), reverse);
  }

  //check for 0 at breakpoint - use when the signal is constant (ie no unique minimum)
  protected void checkBreak(final SamCounts sa, final Distribution di, final boolean reverse, final double tolerance) {
    final Signal sig = new SignalDistributionLn(sa, di, "test");
    final int br = (int) (BREAK + (reverse ? STATS.meanLength() - 1 : 0));
    final double v = sig.value(br);
    //System.err.println("br=" + br + " v=" + v);
    assertEquals(0.0, v, tolerance);
  }

  protected void checkDistr(final SamCounts sa, final Distribution di, final boolean reverse) {
    final Signal sig = new SignalDistributionLn(sa, di, "test");
    checkSignal(sig, reverse);
  }

  private void checkSignal(final Signal sig, final boolean reverse) {
    int mini = -1;
    double minv = Double.POSITIVE_INFINITY;
    for (int i = 0; i < 2 * BREAK; i++) {
      final double v = sig.value(i);
      //System.err.println(i + "\t" + com.rtg.util.Utils.realFormat(v, 2));
      if (v < minv) {
        minv = v;
        mini = i;
      }
    }

    final int diff = (int) (BREAK + (reverse ? STATS.meanLength() : 0)) - mini;
    final String str = "break=" + BREAK + " mini=" + mini + " diff=" + diff + " BREAK_TOLERANCE=" + BREAK_TOLERANCE + " minv=" + minv;
    //System.err.println(str);
    assertTrue(str, -BREAK_TOLERANCE <= diff && diff <= BREAK_TOLERANCE);
  }

  static Signal plotSignal(final AllCounts acs, final BayesianSignal bayesSig, final boolean reverse) {
    final ReadGroupState rgs = new ReadGroupState(STATS, acs);
    return bayesSig.makeSignal(new ReadGroupState[] {rgs}, reverse, "test");
  }

  protected abstract Signal localSignal(final boolean reverse, final boolean debug);

  static void plotSignals(final PrintStream ps, final String label, final AllCounts acs) {
    ps.println("#" + label);
    ps.println("#" + TAB + "delete-left" + TAB + "delete-right" + TAB + "dup-left" + TAB + "dup-right" + TAB + "break" + TAB + "novel" + TAB + "delete" + TAB + "normal" + TAB + "double");
    final Signal[] sigs = {
        plotSignal(acs, new DeleteBoundaryBayesianSignal(), false),
        plotSignal(acs, new DeleteBoundaryBayesianSignal(), true),
        plotSignal(acs, new DuplicateDonorBayesianSignal(), false),
        plotSignal(acs, new DuplicateDonorBayesianSignal(), true),
        plotSignal(acs, new BreakpointBayesianSignal(), false),
        plotSignal(acs, new NovelInsertionBayesianSignal(), false),
        plotSignal(acs, new NormalBayesianSignal(0), false),
        plotSignal(acs, new NormalBayesianSignal(2), false),
        plotSignal(acs, new NormalBayesianSignal(4), false),
    };

    for (int i = -WINDOW_LO; i < BREAK + WINDOW; i++) {
      ps.print(i);
      for (Signal sig : sigs) {
        final double v = sig.value(i);
        ps.print(TAB + Utils.realFormat(v, 2));
      }
      ps.println();
    }
  }

}
