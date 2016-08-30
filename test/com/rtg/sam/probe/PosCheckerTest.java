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
}