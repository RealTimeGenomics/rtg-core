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
package com.rtg.variant;

import java.util.Arrays;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMSequenceRecord;

import junit.framework.TestCase;

/**
 */
public class VariantAlignmentRecordTest extends TestCase {

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
    final VariantAlignmentRecord r = new VariantAlignmentRecord(rec);
    assertTrue(Arrays.equals("TATT".getBytes(), r.getRead()));
    assertEquals(41, r.getStart());
    assertEquals(6, r.getLength());
    assertTrue(Arrays.equals(q, r.getQuality()));
    assertEquals("2M2D2M", r.getCigar());
    assertEquals(43, r.getMappingQuality());
    assertEquals(2, r.getNHOrIH());
    assertEquals("41 2M2D2M TATT ABCD", r.toString());
    assertEquals(false, r.isNegativeStrand());
    assertNull(r.chain());
    r.setNextInChain(r);
    assertEquals(r, r.chain());
    assertFalse(r.isOverflow());
    assertTrue(VariantAlignmentRecord.overflow(0, 1).isOverflow());
    assertFalse(r.isUnmapped());
    assertEquals(-232, r.getFragmentLength());
  }

  public void testNegStrand() {
    final SAMFileHeader header = new SAMFileHeader();
    final SAMRecord rec = new SAMRecord(header);
    rec.setReadNegativeStrandFlag(true);
    VariantAlignmentRecord r = new VariantAlignmentRecord(rec);
    assertTrue(r.isNegativeStrand());
    rec.setReadNegativeStrandFlag(false);
    r = new VariantAlignmentRecord(rec);
    assertFalse(r.isNegativeStrand());
    rec.setReadPairedFlag(true);
    r = new VariantAlignmentRecord(rec);
    assertFalse(r.isNegativeStrand());
    rec.setReadNegativeStrandFlag(true);
    r = new VariantAlignmentRecord(rec);
    assertTrue(r.isNegativeStrand());
    assertFalse(r.isUnmapped());
  }

  public void testUnmapped() {
    final SAMFileHeader header = new SAMFileHeader();
    final SAMRecord rec = new SAMRecord(header);
    rec.setReadUnmappedFlag(true);
    VariantAlignmentRecord r = new VariantAlignmentRecord(rec);
    assertTrue(r.isUnmapped());
  }

  private void check(final VariantAlignmentRecord a, final VariantAlignmentRecord b) {
    assertTrue(a.compareTo(b) > 0);
    assertTrue(b.compareTo(a) < 0);
    assertEquals(0, a.compareTo(a));
  }

  public void testComparator() {
    final SAMFileHeader header = new SAMFileHeader();
    header.addSequence(new SAMSequenceRecord("seq0", 1000));
    header.addSequence(new SAMSequenceRecord("seq1", 1000));
    header.addSequence(new SAMSequenceRecord("seq2", 1000));
    final SAMRecord rec1 = new SAMRecord(header);
    final SAMRecord rec2 = new SAMRecord(header);
    rec1.setReferenceIndex(-1);
    rec2.setReferenceIndex(1);
    check(new VariantAlignmentRecord(rec1), new VariantAlignmentRecord(rec2));
    rec1.setReferenceIndex(2);
    check(new VariantAlignmentRecord(rec1), new VariantAlignmentRecord(rec2));
    rec1.setReferenceIndex(1);
    rec1.setAlignmentStart(11);
    rec2.setAlignmentStart(10);
    check(new VariantAlignmentRecord(rec1), new VariantAlignmentRecord(rec2));
    rec1.setAlignmentStart(10);
    rec1.setCigarString("11M");
    rec2.setCigarString("10M");
    check(new VariantAlignmentRecord(rec1), new VariantAlignmentRecord(rec2));
    rec2.setCigarString("5M6X");
    check(new VariantAlignmentRecord(rec1), new VariantAlignmentRecord(rec2));
    rec2.setCigarString("11M");
    rec1.setMappingQuality(10);
    rec2.setMappingQuality(11);
    check(new VariantAlignmentRecord(rec1), new VariantAlignmentRecord(rec2));
    rec2.setMappingQuality(10);
    // there are more fields for even finer disambiguation
  }

}
