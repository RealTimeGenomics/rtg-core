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

import com.rtg.util.intervals.RegionRestriction;

import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.SAMSequenceRecord;

import junit.framework.TestCase;

/**
 * Test class
 */
public class SamRegionRestrictionTest extends TestCase {

  public void testSomeMethod() {
    SamRegionRestriction srr = new SamRegionRestriction("fooo");
    SAMSequenceDictionary sdd = new SAMSequenceDictionary();
    sdd.addSequence(new SAMSequenceRecord("fooo", 9876));
    assertEquals(0, srr.resolveStart());
    assertEquals(9876, srr.resolveEnd(sdd));
    assertEquals("fooo", srr.getSequenceName());
    assertEquals(RegionRestriction.MISSING, srr.getStart());
    assertEquals(RegionRestriction.MISSING, srr.getEnd());
    assertEquals("fooo", srr.toString());

    srr = new SamRegionRestriction("fooo", 75, 555);
    assertEquals(75, srr.resolveStart());
    assertEquals(555, srr.resolveEnd(sdd));
    assertEquals("fooo", srr.getSequenceName());
    assertEquals(75, srr.getStart());
    assertEquals(555, srr.getEnd());
    assertEquals("fooo:76-555", srr.toString());
    srr = new SamRegionRestriction("fooo", 75, RegionRestriction.MISSING);
    assertEquals("fooo:76", srr.toString());
  }

  public void testParsingConstructor() {
    try {
      new SamRegionRestriction("fooo:2-1");
      fail();
    } catch (IllegalArgumentException e) {
      assertEquals("Malformed range in restriction: \"fooo:2-1\"", e.getMessage());
    }
    try {
      new SamRegionRestriction("fooo:-2-3");
      fail();
    } catch (IllegalArgumentException e) {
      assertEquals("Malformed range in restriction: \"fooo:-2-3\"", e.getMessage());
    }
    SamRegionRestriction srr = new SamRegionRestriction("fooo:1-2");
    assertEquals("fooo", srr.getSequenceName());
    assertEquals(0, srr.getStart());
    assertEquals(2, srr.getEnd());
    assertEquals("fooo:1-2", srr.toString());
  }

}
