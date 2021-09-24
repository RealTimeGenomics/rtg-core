/*
 * Copyright (c) 2018. Real Time Genomics Limited.
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
import htsjdk.samtools.SAMSequenceRecord;
import junit.framework.TestCase;

/**
 */
public class SamReaderRecordTest extends TestCase {

  public void testBasicOperation() {
    final SAMFileHeader header = new SAMFileHeader();
    final SAMRecord rec = new SAMRecord(header);
    rec.setReadString("TATT");
    rec.setReadPairedFlag(true);
    rec.setAlignmentStart(42);
    rec.setCigarString("2M2D2M");
    rec.setMappingQuality(43);
    rec.setAttribute("NH", 2);
    final byte[] q = {32, 33, 34, 35};
    rec.setBaseQualities(q);
    rec.setInferredInsertSize(-232);
    final SamReaderRecord r = new SamReaderRecord.SamReaderRecordPopulator().populate(rec);
    assertEquals(41, r.getStart());
    assertEquals(6, r.getLength());
    assertNull(r.chain());
    r.setNextInChain(r);
    assertEquals(r, r.chain());
    assertEquals(-232, r.getFragmentLength());
  }

  private void check(final SamReaderRecord a, final SamReaderRecord b) {
    assertTrue(a.compareTo(b) > 0);
    assertTrue(b.compareTo(a) < 0);
    assertEquals(0, a.compareTo(a));
  }

  public void testComparator() {
    SamReaderRecord.SamReaderRecordPopulator p = new SamReaderRecord.SamReaderRecordPopulator();
    final SAMFileHeader header = new SAMFileHeader();
    header.addSequence(new SAMSequenceRecord("seq0", 1000));
    header.addSequence(new SAMSequenceRecord("seq1", 1000));
    header.addSequence(new SAMSequenceRecord("seq2", 1000));
    final SAMRecord rec1 = new SAMRecord(header);
    final SAMRecord rec2 = new SAMRecord(header);
    rec1.setReferenceIndex(-1);
    rec2.setReferenceIndex(1);
    check(p.populate(rec1), p.populate(rec2));
    rec1.setReferenceIndex(2);
    check(p.populate(rec1), p.populate(rec2));
    rec1.setReferenceIndex(1);
    rec1.setAlignmentStart(11);
    rec2.setAlignmentStart(10);
    check(p.populate(rec1), p.populate(rec2));
    rec1.setAlignmentStart(10);
    rec1.setCigarString("11M");
    rec2.setCigarString("10M");
    check(p.populate(rec1), p.populate(rec2));
    rec2.setCigarString("5M6X");
    check(p.populate(rec1), p.populate(rec2));
    rec2.setCigarString("11M");
    rec1.setMappingQuality(10);
    rec2.setMappingQuality(11);
    check(p.populate(rec1), p.populate(rec2));
    rec2.setMappingQuality(10);
    // there are more fields for even finer disambiguation
  }
}
