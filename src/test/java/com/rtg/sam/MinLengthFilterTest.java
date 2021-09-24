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