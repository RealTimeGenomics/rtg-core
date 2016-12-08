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

import com.rtg.util.ChiSquared;
import com.rtg.util.diagnostic.Diagnostic;

/**
 */
public class DeleteBoundaryBayesianSignalTest extends AbstractBayesianSignalTest {

  /**
   * Generate a delete distribution.
   * Makes all the same assumptions we made in constructing the distributions.
   * @return signals for a deletion.
   */
  private static AllCounts distribution() {
    final AllCounts sas = new AllCounts(2 * BREAK);

    final double d1 = FRAGMENT_MEAN + WIDTH * FRAGMENT_STD_DEV - MAX_ALIGNMENT;
    final int p1 = (int) (BREAK - d1);
    final double d2 = FRAGMENT_MEAN - WIDTH * FRAGMENT_STD_DEV + MAX_ALIGNMENT - READ_LENGTH;
    final int p2 = (int) (BREAK - d2);
    final int p3 = BREAK - READ_LENGTH + MAX_ALIGNMENT + 1;
    //System.err.println("d1=" + d1 + " p1=" + p1 + " d2=" + d2 + " p2=" + p2 + " p3=" + p3);
    assert 0 < p1 && p1 < p2 && p2 < p3 && p3 < 2 * BREAK;

    //Left arms
    for (int i = 0; i < p1; ++i) {
      sas.properLeft().increment(i, PROPER_RATE);
    }
    for (int i = p1; i < p2; ++i) {
      final double del0 = BREAK - FRAGMENT_MEAN + MAX_ALIGNMENT - i;
      final double del1 = i - BREAK + FRAGMENT_MEAN - READ_LENGTH + MAX_ALIGNMENT;
      final double pr = PROPER_RATE * ChiSquared.normal(del0 / FRAGMENT_STD_DEV);
      final double di = PROPER_RATE * ChiSquared.normal(del1 / FRAGMENT_STD_DEV);
      final double un = PROPER_RATE - pr - di;
      //System.err.println("i=" + i + " del0=" + del0 + " pr=" + pr + " del1=" + del1 + " di=" + di + " sum=" + (pr + di + un));

      sas.properLeft().increment(i, pr);
      sas.discordantLeft().increment(i, di);
      sas.unmatedLeft().increment(i, un);
    }
    for (int i = p2; i < p3; ++i) {
      sas.discordantLeft().increment(i, PROPER_RATE);
    }
    /*
    for (int i = p3; i < 2 * BREAK; ++i) {
      //nothing - unmapped;
    }
    */
    //Right arms
    for (int i = 0; i < p3; ++i) {
      sas.properRight().increment(i, PROPER_RATE);
    }
    //discordant and unmated are random rates

    //add the random rates both arms
    for (int i = 0; i < 2 * BREAK; ++i) {
      sas.properLeft().increment(i, PROPER_RANDOM_RATE);
      sas.properRight().increment(i, PROPER_RANDOM_RATE);
      sas.discordantLeft().increment(i, DISCORDANT_RATE);
      sas.discordantRight().increment(i, DISCORDANT_RATE);
      sas.unmatedLeft().increment(i, UNMATED_RATE);
      sas.unmatedRight().increment(i, UNMATED_RATE);
    }
    //sas.plot(System.err);
    return sas;
  }

  public void testLeftArmProper() {
    for (final boolean reverse : new boolean[] {false, true}) {
      final Distribution d = new DeleteBoundaryBayesianSignal().leftArmProper(STATS, reverse);
      //System.err.println(d.dump());
      assertEquals(PROPER_RATE + PROPER_RANDOM_RATE, d.get(rev(WINDOW_LO, reverse)), 0.01);
      final int d1 = (int) (FRAGMENT_MEAN - MAX_ALIGNMENT);
      final int rev = rev(-d1, reverse);
      assertEquals("" + rev, PROPER_RATE * 0.5 + PROPER_RANDOM_RATE, d.get(rev), 0.01);
      assertEquals(PROPER_RANDOM_RATE, d.get(rev(WINDOW_HI - 1, reverse)), 0.01);

      final SamCounts sa = reverse ? distribution().reverse().properRight() : distribution().properLeft();
      checkDistr(sa, d, reverse);
      checkBreak(sa, d, reverse, 0.02);
    }
  }

  public void testLeftArmDiscor() {
    for (final boolean reverse : new boolean[] {false, true}) {
      final Distribution d = new DeleteBoundaryBayesianSignal().leftArmDiscordant(STATS, reverse);
      //System.err.println(d.dump());
      assertEquals(PROPER_RANDOM_RATE, d.get(rev(WINDOW_LO, reverse)), 0.01);
      final int d1 = (int) (FRAGMENT_MEAN - READ_LENGTH + MAX_ALIGNMENT) - 1;
      final int r1 = rev(-d1, reverse);
      assertEquals(PROPER_RATE * 0.5 + PROPER_RANDOM_RATE, d.get(r1), 0.01);
      final int r2 = rev(-(READ_LENGTH - MAX_ALIGNMENT), reverse);
      assertEquals(PROPER_RATE + PROPER_RANDOM_RATE, d.get(r2), 0.01);
      final int r3 = rev(-(READ_LENGTH - MAX_ALIGNMENT - 1), reverse);
      assertEquals(PROPER_RANDOM_RATE, d.get(r3), 0.01);
      assertEquals(PROPER_RANDOM_RATE, d.get(rev(WINDOW_HI - 1, reverse)), 0.01);

      final SamCounts sa = reverse ? distribution().reverse().discordantRight() : distribution().discordantLeft();
      checkDistr(sa, d, reverse);
      checkBreak(sa, d, reverse, 0.05);
    }
  }

  public void testLeftArmUnmated() {
    for (final boolean reverse : new boolean[] {false, true}) {
      final Distribution d = new DeleteBoundaryBayesianSignal().leftArmUnmated(STATS, reverse);
      //System.err.println(d.dump());

      assertEquals(0.001, d.get(rev(WINDOW_LO, reverse)), 0.001);
      final int p = DistributionTestUtils.peak(d);
      final int r1 = rev((int) -(FRAGMENT_MEAN - READ_LENGTH / 2) + (reverse ? 1 : 0), reverse);
      assertEquals(r1, p);
      assertEquals(0.4182, d.get(p), 0.001); //regression
      assertEquals(0.001, d.get(rev(WINDOW_HI - 1, reverse)), 0.001);

      final SamCounts sa = reverse ? distribution().reverse().unmatedRight() : distribution().unmatedLeft();
      checkDistr(sa, d, reverse);
      checkBreak(sa, d, reverse, 0.07);
    }
  }

  public void testRightArmProper() {
    for (final boolean reverse : new boolean[] {false, true}) {
      final Distribution d = new DeleteBoundaryBayesianSignal().rightArmProper(STATS, reverse);
      //System.err.println(d.dump());
      assertEquals(PROPER_RATE + PROPER_RANDOM_RATE, d.get(rev(WINDOW_LO, reverse)), 0.01);
      final int d1 = -(READ_LENGTH - MAX_ALIGNMENT);
      assertEquals(PROPER_RATE + PROPER_RANDOM_RATE, d.get(rev(d1, reverse)), 0.01);
      assertEquals(PROPER_RANDOM_RATE, d.get(rev(d1 + 1, reverse)), 0.01);
      assertEquals(PROPER_RANDOM_RATE, d.get(rev(WINDOW_HI - 1, reverse)), 0.01);

      final SamArray sa = reverse ? distribution().reverse().properLeft() : distribution().properRight();
      checkDistr(sa, d, reverse);
      checkBreak(sa, d, reverse, 0.01);
    }
  }

  public void testRightArmDiscor() {
    for (final boolean reverse : new boolean[] {false, true}) {
      final Distribution d = new DeleteBoundaryBayesianSignal().rightArmDiscordant(STATS, reverse);
      //System.err.println(d.dump());
      assertEquals(DISCORDANT_RATE, d.get(WINDOW_LO), 0.01);
      assertEquals(DISCORDANT_RATE, d.get(0), 0.01);
      assertEquals(DISCORDANT_RATE, d.get(WINDOW_HI - 1), 0.01);

      final SamCounts sa = reverse ? distribution().reverse().discordantLeft() : distribution().discordantRight();
      checkBreak(sa, d, reverse, 0.01);
    }
  }

  public void testRightArmUnmated() {
    for (final boolean reverse : new boolean[] {false, true}) {
      final Distribution d = new DeleteBoundaryBayesianSignal().rightArmUnmated(STATS, reverse);
      //System.err.println(d.dump());
      assertEquals(UNMATED_RATE, d.get(WINDOW_LO), 0.01);
      assertEquals(UNMATED_RATE, d.get(0), 0.01);
      assertEquals(UNMATED_RATE, d.get(WINDOW_HI - 1), 0.01);

      final SamArray sa = reverse ? distribution().reverse().unmatedLeft() : distribution().unmatedRight();
      checkBreak(sa, d, reverse, 0.01);
    }
  }

  @Override
  protected Signal localSignal(final boolean reverse, final boolean debug) {
    final AllCounts acs = distribution();
    final ReadGroupState rgs = new ReadGroupState(STATS, reverse ? acs.reverse() : acs);
    return new DeleteBoundaryBayesianSignal(debug).makeSignal(new ReadGroupState[] {rgs}, reverse, "test");
  }

  /**
   * Plot the response of all the non-trivial Bayesian signals to this case.
   * @param args ignored
   */
  public static void main(String[] args) {
    Diagnostic.setLogStream();
    plotSignals(System.out, "delete", distribution());
  }

}
