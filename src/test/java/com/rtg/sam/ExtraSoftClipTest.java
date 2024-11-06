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

package com.rtg.sam;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import junit.framework.TestCase;

public class ExtraSoftClipTest extends TestCase {


  public void testExtraSoftClip() {
    check("AAAAAAAAAAAACGTCGATG", "20=", "20=", false, 10);
    check("AAAAAAAAAAAACGTCGATG", "5=1X10=1X3=", "5=1X10=1X3=", false, 10);
    check("AAAAAAAAAAAACGTCGATG", "5=1X10=3X1=", "5=1X10=3X1=", false, 10);
    check("AAAAAAAAAAAACGTCGATG", "2=1X13=3X1=", "3S13=3X1=", false, 13);
    check("AAAAAAAAAAAACGTCGATG", "1X15=3X1=", "1S15=3X1=", false, 11);
    check("AAAAAAAAAAAACGTCGATG", "2=4I10=1X3=", "6S10=1X3=", false, 16);
    check("AAAAAAAAAAAACGTCGATG", "2=4D14=1X3=", "2S14=1X3=", false, 12);
    check("AAAAAAAAAAAACGTCGATG", "20=", "20=", true, 10);
    check("AAAAAAAAAAAACGTCGATG", "5=1X10=1X3=", "5=1X10=1X3=", true, 10);
    check("AAAAAAAAAAAACGTCGATG", "5=1X10=2X2=", "5=1X10=4S", true, 10);
    check("AAAAAAAAAAAACGTCGATG", "5=1X12=2X", "5=1X12=2S", true, 10);

    check("AAAAAAAAAAAACGTCGATG", "15X1=4X", "*", true, 10);
    check("AAAAAAAAAAAACGTCGATG", "15X1=4X", "*", false, 10);
  }

  private void check(String readString, String cigarString, String expReadString, boolean pos, int newStart) {
    final SAMRecord samRecord = new SAMRecord(new SAMFileHeader());
    samRecord.setReadString(readString);
    samRecord.setBaseQualityString(readString);
    samRecord.setCigarString(cigarString);
    samRecord.setReadNegativeStrandFlag(!pos);
    samRecord.setAlignmentStart(10);
    ExtraSoftClip.addSoftClip(samRecord);
    assertEquals(expReadString, samRecord.getCigarString());
    if (!samRecord.getReadUnmappedFlag()) {
      assertEquals(readString.length(), samRecord.getCigar().getReadLength());
    }
    assertEquals(newStart, samRecord.getAlignmentStart());
  }

}
