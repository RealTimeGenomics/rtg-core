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
package com.rtg.util.intervals;

import junit.framework.TestCase;

/**
 */
public class RegionRestrictionTest extends TestCase {

  public void test() {
    checkMalformed("");
    checkMalformed("blah:");
    checkMalformed("blah:600-400");
    checkMalformed("blah:600-+400");
    checkMalformed("blah:600-599");
    checkMalformed("blah:600+0");
    checkMalformed("blah:0+10");
    checkMalformed("chr16:2194033+5000QRST");
    checkMalformed("chr16:2194033-2194060QRST");
    checkMalformed("chr16:2Q-3T");
    checkMalformed("chr16:2-T3");
    checkSuccess("blah:600-600", "blah", 599, 600, "blah:600-600");
    checkSuccess("blah", "blah", RegionRestriction.MISSING, RegionRestriction.MISSING, "blah");
    checkSuccess("blah:1-2", "blah", 0, 2, "blah:1-2");
    checkSuccess("blah:1+2", "blah", 0, 2, "blah:1-2");
    checkSuccess("blah:1+1", "blah", 0, 1, "blah:1-1");
    checkSuccess("blah:1", "blah", 0, -1, "blah:1");

    // These tests assume the default locale uses ',' for thousands and '.' for decimal, make better if need be
    checkMalformed("blah:1000.1-2000");
    checkMalformed("blah:1000-2000.8");
    checkSuccess("blah:1,000-2,000", "blah", 999, 2000, "blah:1000-2000");
  }

  public void testBoundedOverlaps() {
    checkBoundedOverlaps(new RegionRestriction("blah:10-100"));
  }

  // This method currently checks both RegionRestriction and SequenceNameLocusSimple implementations
  static void checkBoundedOverlaps(SequenceNameLocus r) {
    assertEquals(9, r.getStart());
    assertEquals(100, r.getEnd());

    assertFalse(r.overlaps(new SequenceNameLocusSimple("blah", 0, 9)));

    assertTrue(r.overlaps(new SequenceNameLocusSimple("blah", 0, 10)));
    assertTrue(r.overlaps(new SequenceNameLocusSimple("blah", 0, 100)));
    assertTrue(r.overlaps(new SequenceNameLocusSimple("blah", 0, 1000)));
    assertTrue(r.overlaps(new SequenceNameLocusSimple("blah", 11, 20)));
    assertTrue(r.overlaps(new SequenceNameLocusSimple("blah", 9, 1000)));
    assertTrue(r.overlaps(new SequenceNameLocusSimple("blah", 99, 1000)));

    //assertTrue(r.overlaps(new SequenceNameLocusSimple("blah", 9, 9)));
    //assertTrue(r.overlaps(new SequenceNameLocusSimple("blah", 99, 99)));

    assertFalse(r.overlaps(new SequenceNameLocusSimple("blah1", 99, 1000)));
    assertFalse(r.overlaps(new SequenceNameLocusSimple("blah", 100, 1000)));
  }

  // Test overlapping of unbounded regions
  public void testUnboundedOverlaps() {
    RegionRestriction r = new RegionRestriction("blah:10");
    assertFalse(r.overlaps(new SequenceNameLocusSimple("blah", 0, 9)));
    assertTrue(r.overlaps(new SequenceNameLocusSimple("blah", 0, 10)));
    assertTrue(r.overlaps(new SequenceNameLocusSimple("blah", 99, 1000)));

    r = new RegionRestriction("blah");
    assertTrue(r.overlaps(new SequenceNameLocusSimple("blah", 0, 10)));
    assertFalse(r.overlaps(new SequenceNameLocusSimple("blah1", 0, 10)));
  }

  public void testValidateRegion() {
    assertTrue(RegionRestriction.validateRegion("blah"));
    assertFalse(RegionRestriction.validateRegion("blah:"));
  }

  private void checkSuccess(String region, String template, int start, int end, String formatted) {
    checkBasic(template, start, end, formatted);
    final RegionRestriction parsed = new RegionRestriction(region);
    assertEquals(template, parsed.getSequenceName());
    assertEquals(start, parsed.getStart());
    assertEquals(end, parsed.getEnd());
    assertEquals(formatted, parsed.toString());
  }

  private void checkBasic(String template, int start, int end, String formatted) {
    final RegionRestriction expected = new RegionRestriction(template, start, end);
    assertEquals(template, expected.getSequenceName());
    assertEquals(start, expected.getStart());
    assertEquals(end, expected.getEnd());
    assertEquals(formatted, expected.toString());
  }

  private void checkMalformed(String region) {
    try {
      new RegionRestriction(region);
      fail("Expected parsing to grumble about: " + region);
    } catch (IllegalArgumentException e) {
      assertEquals("Malformed range in restriction: \"" + region + "\"", e.getMessage());
    }
  }
}
