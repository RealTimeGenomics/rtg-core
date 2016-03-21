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
package com.rtg.variant;

import com.rtg.ngs.Arm;

/**
 * Interface for phred score scaler objects.
 */
public interface PhredScaler {

  /**
   * Get a phred score from a binary quality optionally
   * correcting it.
   * @param qual original quality value.
   * @param readPosition position on read of <code>qualChar</code>
   * @param arm For paired end reads, which arm this is. Use {@code Arm.LEFT} if single end.
   * @return the possibly corrected phred score.
   */
  int getScaledPhred(byte qual, int readPosition, Arm arm);
}

