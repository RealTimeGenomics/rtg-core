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
 * Looks at collections of SAM records and computes a numeric
 * signal.
 */
public interface Signal {

  /**
   * Get the value of the signal at the specified position.
   * @param position base position of the signal (0 based).
   * @return the value of the signal.
   */
  double value(int position);

  /**
   * Get the column label for the output files.
   * @return the column label for this signal.
   */
  String columnLabel();
}
