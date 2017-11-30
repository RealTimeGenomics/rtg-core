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
package com.rtg.util;

import com.rtg.launcher.globals.GlobalFlags;
import com.rtg.launcher.globals.CoreGlobalFlags;

/**
 */
public final class MaxShiftUtils {

  private MaxShiftUtils() { }

  /**
   * Default base rate for inserted/deleted bases. This is about double the
   * usual sequencer error rate, to allow for some reads being one or two
   * standard deviations above the average.
   */
  private static final double DEFAULT_INDEL_RATE = 0.01;
  /** default total length of all inserts/deletes allowed in reasonably short reads */
  private static final Integer DEFAULT_INDEL_LENGTH = GlobalFlags.getIntegerValue(CoreGlobalFlags.DEFAULT_INDEL_LENGTH_FLAG);

  /**
   * Compute the maximum distance the start and end of the read are allowed to move along the template after
   * alignment, compared to original hit location.
   * This estimate is based on sequencer error rates, and is insufficent for finding longer indels in reads.
   * Should possibly be the maximum number of indels findable based on max score and alignment penalties - although
   * this may have performance implications.
   * @param rlen maximum read length
   * @return the maximum distance the start and end of the read are allowed to
   *         move along the template .
   */
  public static int calculateDefaultMaxShift(int rlen) {
    //return new MaxShiftFactor(factor)
    return (int) (rlen * DEFAULT_INDEL_RATE + DEFAULT_INDEL_LENGTH);
  }
}
