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

import htsjdk.samtools.SAMRecord;

import junit.framework.TestCase;

/**
 */
public class DuplicateSamFilterTest extends TestCase {

  public void testSingleEnd() {
    final DuplicateSamFilter dsf = new DuplicateSamFilter();

    assertFalse(dsf.acceptRecord(null));

    final SAMRecord rec1 = new SAMRecord(null);

    assertFalse(dsf.acceptRecord(rec1));

    rec1.setReadName("IronSheik");

    assertTrue(dsf.acceptRecord(rec1));

    rec1.setNotPrimaryAlignmentFlag(true);

    final SAMRecord rec2 = new SAMRecord(null);
    rec2.setReadName("Sabu");
    rec2.setNotPrimaryAlignmentFlag(true);

    final SAMRecord rec3 = new SAMRecord(null);
    rec3.setReadName("IronSheik");
    rec3.setNotPrimaryAlignmentFlag(true);

    assertTrue(dsf.acceptRecord(rec1));
    assertTrue(dsf.acceptRecord(rec2));
    assertFalse(dsf.acceptRecord(rec1));
    assertFalse(dsf.acceptRecord(rec3));
  }

  public void testPairedEnd() {
    final DuplicateSamFilter dsf = new DuplicateSamFilter();

    assertFalse(dsf.acceptRecord(null));

    final SAMRecord rec1f = new SAMRecord(null);

    assertFalse(dsf.acceptRecord(rec1f));

    final String firstReadName = "IronSheik";
    final String secondReadName = "Sabu";

    rec1f.setReadName(firstReadName);
    rec1f.setReadPairedFlag(true);
    rec1f.setFirstOfPairFlag(true);

    assertTrue(dsf.acceptRecord(rec1f));

    rec1f.setNotPrimaryAlignmentFlag(true);

    final SAMRecord rec2f = new SAMRecord(null);
    rec2f.setReadPairedFlag(true);
    rec2f.setFirstOfPairFlag(true);
    rec2f.setReadName(secondReadName);
    rec2f.setNotPrimaryAlignmentFlag(true);

    final SAMRecord rec3fdup = new SAMRecord(null);
    rec3fdup.setReadPairedFlag(true);
    rec3fdup.setFirstOfPairFlag(true);
    rec3fdup.setReadName(firstReadName);
    rec3fdup.setNotPrimaryAlignmentFlag(true);

    final SAMRecord rec1s = new SAMRecord(null);
    rec1s.setReadName(firstReadName);
    rec1s.setReadPairedFlag(true);
    rec1s.setFirstOfPairFlag(false);
    rec1s.setNotPrimaryAlignmentFlag(true);

    final SAMRecord rec2s = new SAMRecord(null);
    rec2s.setReadPairedFlag(true);
    rec2s.setFirstOfPairFlag(false);
    rec2s.setReadName(secondReadName);
    rec2s.setNotPrimaryAlignmentFlag(true);

    final SAMRecord rec3sdup = new SAMRecord(null);
    rec3sdup.setReadPairedFlag(true);
    rec3sdup.setFirstOfPairFlag(false);
    rec3sdup.setReadName(firstReadName);
    rec3sdup.setNotPrimaryAlignmentFlag(true);


    assertTrue(dsf.acceptRecord(rec1f));
    assertTrue(dsf.acceptRecord(rec1s));
    assertTrue(dsf.acceptRecord(rec2f));
    assertTrue(dsf.acceptRecord(rec2s));
    assertFalse(dsf.acceptRecord(rec1f));
    assertFalse(dsf.acceptRecord(rec1s));
    assertFalse(dsf.acceptRecord(rec3fdup));
    assertFalse(dsf.acceptRecord(rec3sdup));
  }
}
