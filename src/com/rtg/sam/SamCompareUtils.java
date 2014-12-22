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
package com.rtg.sam;

import net.sf.samtools.SAMRecord;

/**
 * Keep the comparison stuff in one place
 */
public final class SamCompareUtils {

  private SamCompareUtils() { }

  /**
   * Compares SAM records on reference index, start position and if paired, mated vs unmated.
   * @param a first record
   * @param b second record
   * @return -1 if first record comes before second record, 1 if the converse. 0 for equal on compared values
   */
  public static int compareSamRecords(SAMRecord a, SAMRecord b) {
    final int thisRef = a.getReferenceIndex();
    final int thatRef = b.getReferenceIndex();
    if (thisRef == -1 && thatRef != -1) {
      return +1;
    } else if (thatRef == -1 && thisRef != -1) {
      return -1;
    }
    //System.err.println("ref this=" + thisRef + " that=" + thatRef);
    if (thisRef < thatRef) {
      return -1;
    }
    if (thisRef > thatRef) {
      return +1;
    }

    final int thisStart = a.getAlignmentStart();
    final int thatStart = b.getAlignmentStart();
    //System.err.println("start this=" + thisStart + " that=" + thatStart);
    if (thisStart < thatStart) {
      return -1;
    }
    if (thisStart > thatStart) {
      return +1;
    }
    if (a.getReadPairedFlag() && !b.getReadPairedFlag()) {
      return -1;
    } else if (b.getReadPairedFlag() && !a.getReadPairedFlag()) {
      return 1;
    }
    if (a.getReadPairedFlag() && b.getReadPairedFlag()) {
      final boolean thisMated = a.getProperPairFlag();
      final boolean thatMated = b.getProperPairFlag();
      if (thisMated && !thatMated) {
        return -1;
      }
      if (!thisMated && thatMated) {
        return 1;
      }
    }
    final boolean thisUnmapped = a.getReadUnmappedFlag();
    final boolean thatUnmapped = b.getReadUnmappedFlag();
    if (thisUnmapped && !thatUnmapped) {
      return 1;
    }
    if (!thisUnmapped && thatUnmapped) {
      return -1;
    }

    return 0;
  }
}
