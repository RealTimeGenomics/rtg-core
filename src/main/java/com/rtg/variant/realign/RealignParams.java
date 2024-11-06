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

package com.rtg.variant.realign;

import com.rtg.util.machine.MachineType;
import com.rtg.variant.util.arithmetic.PossibilityArithmetic;

/**
 * Priors for the re-alignment classes.
 * These correspond to the machine errors that are expected in the reads.
 *
 */
public interface RealignParams {

  /**
   * Get natural log of probability of starting an insert.
   * @return the log of the probability.
   */
  double insertOpenLn();

  /**
   * Get natural log of probability of extending an insert.
   * @return the log of the probability.
   */
  double insertExtendLn();

  /**
   * Get natural log of probability of starting a delete.
   * @return the log of the probability.
   */
  double deleteOpenLn();

  /**
   * Get natural log of probability of extending a delete.
   * @return the log of the probability.
   */
  double deleteExtendLn();

  /**
   * Get the probability of a match.
   * It is returned as the natural log of the probability.
   * @return probability of a match.
   */
  double matchLn();

  /**
   * Get the probability of a mismatch.
   * It is returned as the natural log of the probability.
   * @return probability of a mismatch.
   */
  double misMatchLn();

  /**
   * @return the machine type for reads, or null if the parameters do not correspond to reads.
   */
  MachineType machineType();

  /**
   * The minimum (most negative) width of a given gap.
   *
   * @param gap the gaps are numbered from 0 (the overlap) to 2 (the large gap).
   * @return For example: -3 for the overlap region (gap == 0).
   * @throws UnsupportedOperationException if <code>completeGenomics()</code> is false.
   */
  int gapStart(int gap);

  /**
   * The maximum (most positive) width of a given gap.
   *
   * @param gap the gaps are numbered from 0 (the overlap) to 2 (the large gap).
   * @return For example: -1 for the overlap region (gap == 0).
   * @throws UnsupportedOperationException if <code>completeGenomics()</code> is false.
   */
  int gapEnd(int gap);

  /**
   * Get the probability of a gap of a given width.
   * It is returned as the natural log of the probability.
   *
   * @param gap the gaps are numbered from 0 (the overlap) to 2 (the large gap).
   * @param width must be between <code>gapStart(gap) .. gapEnd(gap)</code>, inclusive.
   * @return For example: <code>ln(0.84)</code> for a width of -2, when gap == 0.
   * @throws UnsupportedOperationException if <code>completeGenomics()</code> is false.
   */
  double gapFreqLn(int gap, int width);

  /**
   * @param arith arithmetic delegate.
   * @return an array of possibilities (as determined by arith) for each CG gap.
   */
  double[][] gapDistributionPoss(final PossibilityArithmetic arith);
}
