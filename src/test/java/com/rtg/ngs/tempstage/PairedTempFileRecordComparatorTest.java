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
package com.rtg.ngs.tempstage;

import com.rtg.sam.SamBamConstants;

import junit.framework.TestCase;

/**
 */
public class PairedTempFileRecordComparatorTest extends TestCase {

  public void testComparatorDirectly() {
        final PairedTempFileRecordComparator comp = new PairedTempFileRecordComparator();
        final BinaryTempFileRecord samrec1 = makeRecord(5, 1, false, 6, true, true);
        final BinaryTempFileRecord samrec2 = makeRecord(5, 1, false, 6, true, false);
        final BinaryTempFileRecord samrec3 = makeRecord(5, 1, false, 6, true, true);
        final BinaryTempFileRecord samrec4 = makeRecord(5, 1, true, 6, true, true);
        final BinaryTempFileRecord samrec5 = makeRecord(5, 1, true, 6, false, true);
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
                                        final boolean negStrand,
                       final int mateAlignStart, final boolean mateNegStrand, final boolean first) {
    final BinaryTempFileRecord rec = new BinaryTempFileRecord(true, false, false, false);
    rec.setStartPosition(alignStart);
    rec.setReadId(readId);
    rec.setMatePosition(mateAlignStart);
    int samFlags = 0;
    if (negStrand) {
      samFlags |= SamBamConstants.SAM_READ_IS_REVERSE;
    }
      samFlags |= SamBamConstants.SAM_READ_IS_PAIRED;
      if (first) {
        samFlags |= SamBamConstants.SAM_READ_IS_FIRST_IN_PAIR;
      } else {
        samFlags |= SamBamConstants.SAM_READ_IS_SECOND_IN_PAIR;
      }
    if (mateNegStrand) {
      samFlags |= SamBamConstants.SAM_MATE_IS_REVERSE;
    }
    rec.setSamFlags((byte) samFlags);
    return rec;
  }
}
