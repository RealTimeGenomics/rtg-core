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

package com.rtg.sam.probe;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import junit.framework.TestCase;

/**
 *
 */
public class PosCheckerTest extends TestCase {

  static SAMRecord createRecord(String read, String cigar) {
    final SAMRecord rec = new SAMRecord(new SAMFileHeader());
    rec.setReferenceName("A");
    rec.setReadName("read");
    rec.setReadString(read);
    rec.setBaseQualityString(String.format("%" + read.length() + "s", "").replace(' ', '#'));
    rec.setCigarString(cigar);
    rec.setAlignmentStart(1001);
    return rec;
  }

  public void testSetAlignment() {
    final PosChecker pos = new PosChecker(10);
    final SAMRecord rec = createRecord("AGGTTTGG", "1=1X1=1X1=1X1=1X");
    pos.setAlignmentStart(rec, null, 1002);
    assertEquals("GTTTGG", rec.getReadString());
    assertEquals("1=1X1=1X1=1X", rec.getCigarString());
  }

  public void testTrimAllProbe() {
    final PosChecker pos = new PosChecker(10);
    final SAMRecord rec = createRecord("AGGTTTGG", "1=1X1=1X1=1X1=1X");
    pos.setAlignmentStart(rec, null, 2002);
    assertEquals("*", rec.getReadString());
    assertEquals("*", rec.getCigarString());
  }
}
