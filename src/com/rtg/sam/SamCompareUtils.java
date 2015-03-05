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

import htsjdk.samtools.SAMRecord;

/**
 * Keep the comparison stuff in one place
 */
public final class SamCompareUtils {

  private SamCompareUtils() { }

  private static int reverse3bits(final int x) {
    // Take last three bits abc of x and return bca
    return (0b111_011_101_001_110_010_100_000 >>> ((x & 7) * 3)) & 7;
  }

  /**
   * Compares SAM records on reference index, start position and if paired, mated vs unmated.
   * @param a first record
   * @param b second record
   * @return -1 if first record comes before second record, 1 if the converse. 0 for equal on compared values
   */
  public static int compareSamRecords(SAMRecord a, SAMRecord b) {
    final int thisRef = a.getReferenceIndex() & 0x7FFFFFFF; // makes -1 a bignum
    final int thatRef = b.getReferenceIndex() & 0x7FFFFFFF;
    final int rc = Integer.compare(thisRef, thatRef);
    if (rc != 0) {
      return rc;
    }

    final int ac = Integer.compare(a.getAlignmentStart(), b.getAlignmentStart());
    if (ac != 0) {
      return ac;
    }

    // Do this ... (this one doesn't have same ordering as original)
    //return Integer.compare(~a.getFlags(), ~b.getFlags());

    // or this ...
    // Following compares READ_PAIRED_FLAG, PORPER_PAIRED_FLAG, UNMAPPED_FLAG in that order.
    // The ^3 toggles values to get the ordering we require
    // The reverse makes sure comparison is in the order we require
    return Integer.compare(reverse3bits(a.getFlags() ^ 3), reverse3bits(b.getFlags() ^ 3));

    // or this ...

//    final int rpc = Boolean.compare(b.getReadPairedFlag(), a.getReadPairedFlag());
//    if (rpc != 0) {
//      return rpc;
//    }
//
//    if (a.getReadPairedFlag()) {
//      assert b.getReadPairedFlag();
//      final int mc = Boolean.compare(b.getProperPairFlag(), a.getProperPairFlag());
//      if (mc != 0) {
//        return mc;
//      }
//    }
//
//    return Boolean.compare(a.getReadUnmappedFlag(), b.getReadUnmappedFlag());
  }
}
