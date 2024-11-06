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
