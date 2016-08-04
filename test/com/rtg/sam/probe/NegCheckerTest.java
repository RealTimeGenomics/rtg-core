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

import htsjdk.samtools.SAMRecord;
import junit.framework.TestCase;

/**
 *
 */
public class NegCheckerTest extends TestCase {

  public void test() {
    final NegChecker neg = new NegChecker(10);
    final SAMRecord rec = PosCheckerTest.createRecord("AGGTTTGG", "1=1X1=1X3=1X");
    neg.setAlignmentEnd(rec, 1006);
    assertEquals("AGGTTT", rec.getReadString());
    assertEquals("1=1X1=1X2=", rec.getCigarString());
  }
}