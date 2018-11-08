/*
 * Copyright (c) 2016. Real Time Genomics Limited.
 *
 * Use of this source code is bound by the Real Time Genomics Limited Software Licence Agreement
 * for Academic Non-commercial Research Purposes only.
 *
 * If you did not receive a license accompanying this file, a copy must first be obtained by email
 * from support@realtimegenomics.com.  On downloading, using and/or continuing to use this source
 * code you accept the terms of that license agreement and any amendments to those terms that may
 * be made from time to time by Real Time Genomics Limited.
 */
package com.rtg.sam;

import htsjdk.samtools.SAMRecord;

/**
 * Look for cases that were not soft-clipped during primary alignment but perhaps should have been.
 */
public final class MinLengthFilter {

  private MinLengthFilter() { }

  /**
   * Masks mismatches and inserts present at the end of the strand in given read
   * @param record sam record to modify
   * @param minReadLength length of the shortest read that will be kept
   * @return true if record was modified
   */
  public static boolean filterShortReads(SAMRecord record, int minReadLength) {
    if (record.getAlignmentEnd() - record.getAlignmentStart() + 1 < minReadLength) {
      SamUtils.convertToUnmapped(record);
      return true;
    }
    return false;
  }
}
