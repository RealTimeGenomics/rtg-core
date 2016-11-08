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
    //check("AAAAAAAAAAAACGTCGATG", "15X1=4X", "20S", false); // ??
    check("AAAAAAAAAAAACGTCGATG", "20=", "20=", true, 10);
    check("AAAAAAAAAAAACGTCGATG", "5=1X10=1X3=", "5=1X10=1X3=", true, 10);
    check("AAAAAAAAAAAACGTCGATG", "5=1X10=2X2=", "5=1X10=4S", true, 10);
    check("AAAAAAAAAAAACGTCGATG", "5=1X12=2X", "5=1X12=2S", true, 10);
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
    assertEquals(readString.length(), samRecord.getCigar().getReadLength());
    assertEquals(newStart, samRecord.getAlignmentStart());
  }

}