/*
 * Copyright (c) 2018. Real Time Genomics Limited.
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
