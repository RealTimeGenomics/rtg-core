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
package com.rtg.variant.bayes.multisample;

/**
 * Measure the number of repetitive bases between two points on a template.
 */
public interface RepeatMeasurer {

  /**
   * Measures the total length of simple repeat regions between two positions.
   * @param positionA template start position
   * @param positionB template end position (exclusive)
   * @return the total length of repeats
   */
  int measureRepeats(int positionA, int positionB);

  /**
   * Measures the total length of simple repeat regions between two positions.
   * @param positionA template start position
   * @param positionB template end position (exclusive)
   * @param repeatHint hint as to the length of repeat unit that may be present
   * @return the total length of repeats
   */
  int measureRepeats(int positionA, int positionB, int repeatHint);
}
