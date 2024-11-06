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

import com.rtg.util.ChiSquared;
import com.rtg.util.diagnostic.Diagnostic;

/**
 */
public class DuplicateDonorBayesianSignalTest extends AbstractBayesianSignalTest {

  /**
   * Generate a duplicate donor distribution.
   * Makes all the same assumptions we made in constructing the distributions.
   * @return signals for a deletion.
   */
  static AllCounts distribution() {
    final AllCounts sas = new AllCounts(2 * BREAK);

    final int p1 = BREAK - MAX_ALIGNMENT + 1;
    final double d2 = FRAGMENT_MEAN - WIDTH * FRAGMENT_STD_DEV + MAX_ALIGNMENT - 2 * READ_LENGTH;
    final int p2 = (int) (BREAK + d2);
    final double d3 = FRAGMENT_MEAN + WIDTH * FRAGMENT_STD_DEV - MAX_ALIGNMENT - READ_LENGTH;
    final int p3 = (int) (BREAK + d3);
    assert 0 < p1 && p1 < p2 && p2 < p3 && p3 < 2 * BREAK;

    //Left arms
    for (int i = 0; i < p1; ++i) {
      sas.properLeft().increment(i, PROPER_RATE);
    }
    for (int i = p1; i < 2 * BREAK; ++i) {
      sas.properLeft().increment(i, 2 * PROPER_RATE);
      sas.unmatedLeft().increment(i, UNMATED_RATE);
      sas.discordantLeft().increment(i, DISCORDANT_RATE);
    }

    //Right arms
    for (int i = 0; i < p1; ++i) {
      sas.properRight().increment(i, PROPER_RATE);
    }
    for (int i = p1; i < p2; ++i) {
      sas.properRight().increment(i, PROPER_RATE);
      sas.discordantRight().increment(i, PROPER_RATE);
    }
    for (int i = p2; i < p3; ++i) {
      final double delPr = i - BREAK - (FRAGMENT_MEAN - MAX_ALIGNMENT - READ_LENGTH);
      final double pr = PROPER_RATE * ChiSquared.normal(delPr / FRAGMENT_STD_DEV);
      final double delDi = FRAGMENT_MEAN - 2 * READ_LENGTH + MAX_ALIGNMENT - (i - BREAK);
      final double di = PROPER_RATE * ChiSquared.normal(delDi / FRAGMENT_STD_DEV);
      final double un = PROPER_RATE - pr - di;
      //System.err.println("i=" + i + " del0=" + del0 + " pr=" + pr + " del1=" + del1 + " di=" + di + " sum=" + (pr + di + un));

      sas.properRight().increment(i, pr + PROPER_RATE);
      sas.discordantRight().increment(i, di);
      sas.unmatedRight().increment(i, un);
    }
    for (int i = p3; i < 2 * BREAK; ++i) {
      sas.properRight().increment(i, 2 * PROPER_RATE);
    }

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
      final Distribution d = new DuplicateDonorBayesianSignal().leftArmProper(STATS, reverse);
      //System.err.println(d.dump());
      assertEquals(PROPER_RATE + PROPER_RANDOM_RATE, d.get(rev(WINDOW_LO, reverse)), 0.01);
      assertEquals(PROPER_RATE + PROPER_RANDOM_RATE, d.get(rev(-MAX_ALIGNMENT, reverse)), 0.01);
      assertEquals(2 * PROPER_RATE + PROPER_RANDOM_RATE, d.get(rev(-MAX_ALIGNMENT + 1, reverse)), 0.01);
      assertEquals(2 * PROPER_RATE + PROPER_RANDOM_RATE, d.get(rev(WINDOW_HI - 1, reverse)), 0.01);

      final SamCounts sa = reverse ? distribution().reverse().properRight() : distribution().properLeft();
      checkDistr(sa, d, reverse);
      checkBreak(sa, d, reverse, 0.01);
    }
  }

  public void testLeftArmDiscor() {
    for (final boolean reverse : new boolean[] {false, true}) {
      final Distribution d = new DuplicateDonorBayesianSignal().leftArmDiscordant(STATS, reverse);
      //System.err.println(d.dump());
      assertEquals(DISCORDANT_RATE, d.get(WINDOW_LO), 0.01);
      assertEquals(DISCORDANT_RATE, d.get(0), 0.01);
      assertEquals(DISCORDANT_RATE, d.get(WINDOW_HI - 1), 0.01);

      final SamCounts sa = reverse ? distribution().reverse().discordantRight() : distribution().discordantLeft();
      checkBreak(sa, d, reverse, 0.01);
    }
  }

  public void testLeftArmUnmated() {
    for (final boolean reverse : new boolean[] {false, true}) {
      final Distribution d = new DuplicateDonorBayesianSignal().leftArmUnmated(STATS, reverse);
      //System.err.println(d.dump());
      assertEquals(UNMATED_RATE, d.get(rev(WINDOW_LO, reverse)), 0.01);
      assertEquals(UNMATED_RATE, d.get(rev(-MAX_ALIGNMENT, reverse)), 0.01);
      assertEquals(2 * UNMATED_RATE, d.get(rev(-MAX_ALIGNMENT + 1, reverse)), 0.01);
      assertEquals(2 * UNMATED_RATE, d.get(rev(WINDOW_HI - 1, reverse)), 0.01);

      final SamCounts sa = reverse ? distribution().reverse().unmatedRight() : distribution().unmatedLeft();
      checkBreak(sa, d, reverse, 0.01);
    }
  }

  public void testRightArmProper() {
    for (final boolean reverse : new boolean[] {false, true}) {
      final Distribution d = new DuplicateDonorBayesianSignal().rightArmProper(STATS, reverse);
      //System.err.println(d.dump());
      assertEquals(PROPER_RATE + PROPER_RANDOM_RATE, d.get(rev(WINDOW_LO, reverse)), 0.01);
      final int rev = rev((int) (FRAGMENT_MEAN - READ_LENGTH - MAX_ALIGNMENT) + 1, reverse);
      assertEquals("" + rev, PROPER_RATE * 1.5 + PROPER_RANDOM_RATE, d.get(rev), 0.01);
      assertEquals(2 * PROPER_RATE + PROPER_RANDOM_RATE, d.get(rev(WINDOW_HI - 1, reverse)), 0.01);

      final SamCounts sa = reverse ? distribution().reverse().properLeft() : distribution().properRight();
      checkDistr(sa, d, reverse);
      checkBreak(sa, d, reverse, 0.1);
    }
  }

  public void testRightArmDiscor() {
    for (final boolean reverse : new boolean[] {false, true}) {
      final Distribution d = new DuplicateDonorBayesianSignal().rightArmDiscordant(STATS, reverse);
      //System.err.println(d.dump());
      assertEquals(DISCORDANT_RATE, d.get(rev(WINDOW_LO, reverse)), 0.01);
      assertEquals(DISCORDANT_RATE, d.get(rev(-MAX_ALIGNMENT, reverse)), 0.01);
      assertEquals(PROPER_RATE + PROPER_RANDOM_RATE, d.get(rev(-MAX_ALIGNMENT + 1, reverse)), 0.01);
      final int rev = rev((int) (FRAGMENT_MEAN - 2 * READ_LENGTH + MAX_ALIGNMENT), reverse);
      assertEquals("" + rev, PROPER_RATE * 0.5 + PROPER_RANDOM_RATE, d.get(rev), 0.01);
      assertEquals(DISCORDANT_RATE, d.get(WINDOW_HI - 1), 0.01);

      final SamCounts sa = reverse ? distribution().reverse().discordantLeft() : distribution().discordantRight();
      checkDistr(sa, d, reverse);
      checkBreak(sa, d, reverse, 0.02);
    }
  }

  public void testRightArmUnmated() {
    for (final boolean reverse : new boolean[] {false, true}) {
      final Distribution d = new DuplicateDonorBayesianSignal().rightArmUnmated(STATS, reverse);
      //System.err.println(d.dump());
      assertEquals(UNMATED_RATE, d.get(rev(WINDOW_LO, reverse)), 0.01);
      assertEquals(UNMATED_RATE, d.get(rev(-MAX_ALIGNMENT, reverse)), 0.01);
      final int rev = rev((int) (FRAGMENT_MEAN - 3 * READ_LENGTH / 2), reverse);
      assertEquals("" + rev, 0.4182, d.get(rev), 0.01); //regression
      assertEquals(UNMATED_RATE, d.get(rev(WINDOW_HI - 1, reverse)), 0.01);

      final SamCounts sa = reverse ? distribution().reverse().unmatedLeft() : distribution().unmatedRight();
      checkDistr(sa, d, reverse);
      checkBreak(sa, d, reverse, 0.07);
    }
  }

  @Override
  protected Signal localSignal(final boolean reverse, final boolean debug) {
    final AllCounts acs = distribution();
    final ReadGroupState rgs = new ReadGroupState(STATS, reverse ? acs.reverse() : acs);
    return new DuplicateDonorBayesianSignal(debug).makeSignal(new ReadGroupState[] {rgs}, reverse, "test");
  }

  /**
   * Plot the response of all the non-trivial Bayesian signals to this case.
   * @param args ignored
   */
  public static void main(String[] args) {
    Diagnostic.setLogStream();
    plotSignals(System.out, "duplicate", distribution());
  }

}
