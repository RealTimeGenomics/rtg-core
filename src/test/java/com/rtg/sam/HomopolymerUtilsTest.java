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

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import junit.framework.TestCase;

public class HomopolymerUtilsTest extends TestCase {


  public void testMaskHomopolymer() {
    check("AAAAAAAAAAAACGTCGATG", "20=", "AAAAAAAAAAAACGTCGATG", false);
    check("AAAAAAAAAAAACGTCGATG", "5=1X10=1X3=", "AAAAANAAAAAACGTCGATG", false);
    check("AAAAAAAAAAAACGTCGATG", "20=", "AAAAAAAAAAAACGTCGATG", false);
    check("CGTCGATGAAAAAAAAAAAA", "5=1X10=1X3=", "CGTCGATGAAAAAAAANAAA", true);
    check("AAAAAAAAAAAACGTCGATG", "20=", "AAAAAAAAAAAACGTCGATG", true);
    check("AAAAAAAAAAAAAAAAAAAA", "5=1X10=1X3=", "AAAAANAAAAAAAAAANAAA", false);
    check("AAAAAAAAAAAAAAAAAAAA", "5=1X10=1X3=", "AAAAANAAAAAAAAAANAAA", true);
  }

  private void check(String readString, String cigarString, String expReadString, boolean pos) {
    final SAMRecord samRecord = new SAMRecord(new SAMFileHeader());
    samRecord.setReadString(readString);
    samRecord.setBaseQualityString(readString);
    samRecord.setCigarString(cigarString);
    samRecord.setReadNegativeStrandFlag(!pos);
    HomopolymerUtils.maskHomoPolymer(samRecord);
    assertEquals(expReadString, samRecord.getReadString());
  }

}