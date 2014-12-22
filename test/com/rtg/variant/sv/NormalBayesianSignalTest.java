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

/**
 */
public class NormalBayesianSignalTest extends AbstractBayesianSignalTest {

  public void testLeftArmProper() {
    for (final boolean reverse : new boolean[] {false, true}) {
      final Distribution d = new NormalBayesianSignal(4).leftArmProper(STATS, reverse);
      //System.err.println(d.dump());
      final double rate = 2 * PROPER_RATE + PROPER_RANDOM_RATE;
      assertEquals(rate, d.get(rev(WINDOW_LO, reverse)), 0.01);
      assertEquals(rate, d.get(rev(0, reverse)), 0.01);
      assertEquals(rate, d.get(rev(WINDOW_HI - 1, reverse)), 0.01);

      checkBreak(constant(rate), d, reverse, 0.01);
    }
  }

  public void testLeftArmDiscor() {
    for (final boolean reverse : new boolean[] {false, true}) {
      final Distribution d = new NormalBayesianSignal(4).leftArmDiscordant(STATS, reverse);
      //System.err.println(d.dump());
      assertEquals(DISCORDANT_RATE, d.get(rev(WINDOW_LO, reverse)), 0.01);
      assertEquals(DISCORDANT_RATE, d.get(rev(0, reverse)), 0.01);
      assertEquals(DISCORDANT_RATE, d.get(rev(WINDOW_HI - 1, reverse)), 0.01);

      checkBreak(constant(DISCORDANT_RATE), d, reverse, 0.01);
    }
  }

  public void testLeftArmUnmated() {
    for (final boolean reverse : new boolean[] {false, true}) {
      final Distribution d = new NormalBayesianSignal(4).leftArmUnmated(STATS, reverse);
      //System.err.println(d.dump());

      assertEquals(UNMATED_RATE, d.get(rev(WINDOW_LO, reverse)), 0.01);
      assertEquals(UNMATED_RATE, d.get(rev(0, reverse)), 0.01);
      assertEquals(UNMATED_RATE, d.get(rev(WINDOW_HI - 1, reverse)), 0.01);

      checkBreak(constant(UNMATED_RATE), d, reverse, 0.01);
    }
  }

  public void testRightArmProper() {
    for (final boolean reverse : new boolean[] {false, true}) {
      final Distribution d = new NormalBayesianSignal(4).rightArmProper(STATS, reverse);
      //System.err.println(d.dump());
      final double rate = 2 * PROPER_RATE + PROPER_RANDOM_RATE;
      assertEquals(rate, d.get(rev(WINDOW_LO, reverse)), 0.01);
      assertEquals(rate, d.get(rev(0, reverse)), 0.01);
      assertEquals(rate, d.get(rev(WINDOW_HI - 1, reverse)), 0.01);

      checkBreak(constant(rate), d, reverse, 0.01);
    }
  }

  public void testRightArmDiscor() {
    for (final boolean reverse : new boolean[] {false, true}) {
      final Distribution d = new NormalBayesianSignal(4).rightArmDiscordant(STATS, reverse);
      //System.err.println(d.dump());
      assertEquals(DISCORDANT_RATE, d.get(rev(WINDOW_LO, reverse)), 0.01);
      assertEquals(DISCORDANT_RATE, d.get(rev(0, reverse)), 0.01);
      assertEquals(DISCORDANT_RATE, d.get(rev(WINDOW_HI - 1, reverse)), 0.01);

      checkBreak(constant(DISCORDANT_RATE), d, reverse, 0.01);
    }
  }

  public void testRightArmUnmated() {
    for (final boolean reverse : new boolean[] {false, true}) {
      final Distribution d = new NormalBayesianSignal(4).rightArmUnmated(STATS, reverse);
      //System.err.println(d.dump());
      assertEquals(UNMATED_RATE, d.get(rev(WINDOW_LO, reverse)), 0.01);
      assertEquals(UNMATED_RATE, d.get(rev(0, reverse)), 0.01);
      assertEquals(UNMATED_RATE, d.get(rev(WINDOW_HI - 1, reverse)), 0.01);

      checkBreak(constant(UNMATED_RATE), d, reverse, 0.01);
    }
  }

  @Override
  public void testBreakDetection() {
    //no peak to detect
  }

  @Override
  protected Signal localSignal(final boolean reverse, final boolean debug) {
    throw new UnsupportedOperationException();
  }
}
