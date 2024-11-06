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
