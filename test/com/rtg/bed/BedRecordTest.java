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

package com.rtg.bed;

import junit.framework.TestCase;

/**
 */
public class BedRecordTest extends TestCase {

  public void testBedRecord() {
    BedRecord rec = new BedRecord("chr1", 2, 80, "anno1", "anno2");
    assertEquals("chr1", rec.getSequenceName());
    assertEquals(2, rec.getStart());
    assertEquals(80, rec.getEnd());
    assertEquals(2, rec.getAnnotations().length);
    assertEquals("anno1", rec.getAnnotations()[0]);
    assertEquals("anno2", rec.getAnnotations()[1]);
    assertEquals("chr1\t2\t80\tanno1\tanno2", rec.toString());
    rec = new BedRecord("chr1", 2, 80);
    assertNotNull(rec.getAnnotations());
    assertEquals(0, rec.getAnnotations().length);
  }

  public void testBedParsing() {
    final BedRecord rec = BedRecord.fromString("chr1\t2\t80\tanno1\tanno2");
    assertEquals("chr1", rec.getSequenceName());
    assertEquals(2, rec.getStart());
    assertEquals(80, rec.getEnd());
    assertEquals(2, rec.getAnnotations().length);
    assertEquals("anno1", rec.getAnnotations()[0]);
    assertEquals("anno2", rec.getAnnotations()[1]);
    assertEquals("chr1\t2\t80\tanno1\tanno2", rec.toString());
    try {
      BedRecord.fromString("chr1\tads2\t123\tdasf");
      fail();
    } catch (NumberFormatException e) {
      //Expected
    }

    try {
      BedRecord.fromString("chr1");
      fail();
    } catch (ArrayIndexOutOfBoundsException e) {
      //Expected
    }
  }
}
