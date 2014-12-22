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
package com.rtg.ngs.tempstage;

import com.rtg.sam.SamBamConstants;

import junit.framework.TestCase;

/**
 */
public class TempFileRecordComparatorTest extends TestCase {

  public void testComparatorDirectly() {
        final TempFileRecordComparator comp = new TempFileRecordComparator();
        BinaryTempFileRecord samrec1 = makeRecord(5, 1, false);
        BinaryTempFileRecord samrec2 = makeRecord(6, 1, false);
        BinaryTempFileRecord samrec3 = makeRecord(5, 2, false);
        BinaryTempFileRecord samrec4 = makeRecord(5, 1, true);
        BinaryTempFileRecord samrec5 = makeRecord(5, 2, true);
        assertEquals(-1, comp.compare(samrec1, samrec2));
        assertEquals(1, comp.compare(samrec2, samrec1));
        assertEquals(1, comp.compare(samrec2, samrec3));
        assertEquals(-1, comp.compare(samrec3, samrec2));
        assertEquals(1, comp.compare(samrec3, samrec4));
        assertEquals(-1, comp.compare(samrec4, samrec3));
        assertEquals(-1, comp.compare(samrec4, samrec5));
        assertEquals(1, comp.compare(samrec5, samrec4));
  }

  BinaryTempFileRecord makeRecord(final int alignStart, final int readId,
                                        final boolean negStrand) {
    BinaryTempFileRecord rec = new BinaryTempFileRecord(false, false, false, false);
    rec.setStartPosition(alignStart);
    rec.setReadId(readId);
    int samFlags = 0;
    if (negStrand) {
      samFlags |= SamBamConstants.SAM_READ_IS_REVERSE;
    }
    /*
    if (readPaired) {
      samFlags |= SamBamConstants.SAM_READ_IS_PAIRED;
      if (first) {
        samFlags |= SamBamConstants.SAM_READ_IS_FIRST_IN_PAIR;
      } else {
        samFlags |= SamBamConstants.SAM_READ_IS_SECOND_IN_PAIR;
      }
    }
    if (readPaired) {
      samFlags |= SamBamConstants.SAM_READ_IS_MAPPED_IN_PROPER_PAIR;
    }
    if (mateNegStrand) {
      samFlags |= SamBamConstants.SAM_MATE_IS_REVERSE;
    }
    rec.setSamFlags();
    rec.setReadNegativeStrandFlag(negStrand);
    rec.setMateAlignmentStart(mateAlignStart);
    rec.setMateNegativeStrandFlag(mateNegStrand);
    rec.setReadPairedFlag(readPaired);
    rec.setFirstOfPairFlag(first);
    */
    rec.setSamFlags((byte) samFlags);
    return rec;
  }
}
