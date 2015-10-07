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
  int gapStart(int gap) throws UnsupportedOperationException;

  /**
   * The maximum (most positive) width of a given gap.
   *
   * @param gap the gaps are numbered from 0 (the overlap) to 2 (the large gap).
   * @return For example: -1 for the overlap region (gap == 0).
   * @throws UnsupportedOperationException if <code>completeGenomics()</code> is false.
   */
  int gapEnd(int gap) throws UnsupportedOperationException;

  /**
   * Get the probability of a gap of a given width.
   * It is returned as the natural log of the probability.
   *
   * @param gap the gaps are numbered from 0 (the overlap) to 2 (the large gap).
   * @param width must be between <code>gapStart(gap) .. gapEnd(gap)</code>, inclusive.
   * @return For example: <code>ln(0.84)</code> for a width of -2, when gap == 0.
   * @throws UnsupportedOperationException if <code>completeGenomics()</code> is false.
   */
  double gapFreqLn(int gap, int width) throws UnsupportedOperationException;

  /**
   * @param arith arithmetic delegate.
   * @return an array of possibilities (as determined by arith) for each CG gap.
   */
  double[][] gapDistributionPoss(final PossibilityArithmetic arith);
}
