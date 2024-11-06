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

import static org.junit.Assert.assertEquals;

import org.junit.Test;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;

/**
 */
public class MinLengthFilterTest {
  @Test
  public void testMinLength() {
    check("AAAAAAAAAAAACGTCGATG", "20=", true, false);
    check("AAAA", "4=", false, true);
    check("AAAAA", "5=", false, false);
  }

  private void check(String readString, String cigarString, boolean pos, boolean expectedResult) {
    final SAMRecord samRecord = new SAMRecord(new SAMFileHeader());
    samRecord.setReadString(readString);
    samRecord.setBaseQualityString(readString);
    samRecord.setCigarString(cigarString);
    samRecord.setReadNegativeStrandFlag(!pos);
    samRecord.setAlignmentStart(10);
    assertEquals(expectedResult, MinLengthFilter.filterShortReads(samRecord, 5));
    assertEquals(readString, samRecord.getReadString());
    assertEquals(10, samRecord.getAlignmentStart());
  }

}
