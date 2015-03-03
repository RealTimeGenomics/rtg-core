/*
 * Copyright (c) 2015. Real Time Genomics Limited.
 *
 * Use of this source code is bound by the Real Time Genomics Limited Software Licence Agreement
 * for Academic Non-commercial Research Purposes only.
 *
 * If you did not receive a license accompanying this file, a copy must first be obtained by email
 * from support@realtimegenomics.com.  On downloading, using and/or continuing to use this source
 * code you accept the terms of that license agreement and any amendments to those terms that may
 * be made from time to time by Real Time Genomics Limited.
 */
package com.rtg.metagenomics;

/**
 * Interface for different types of root finding mechanisms.
 * Looks for a zero in the function.
 */
public interface LineSolver {

  /**
   * Finds point where line crosses y at 0.  Assumes this is at a point where x is greater than 0.
   * @param line line the process
   * @param relThreshold termination threshold
   * @return cross over point
   */
  double solveLine(Line line, double relThreshold);

}
