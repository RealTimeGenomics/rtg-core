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
package com.rtg.variant.coverage;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceRecord;

import junit.framework.TestCase;

/**
 */
public class DummyMateInfoReaderRecordTest extends TestCase {
  private static class MyRecord extends AbstractMateInfoReaderRecord<MyRecord> {
    MyRecord(SAMRecord sam, int genome) {
      super(sam, genome);
    }

    @Override
    public int disambiguateDuplicate(MyRecord rec) {
      return 0;
    }

    @Override
    public int compareTo(MyRecord o) {
      return 0;
    }
    @Override
    public boolean equals(Object o) {
      if (o == null) {
        return false;
      }
      if (!(o instanceof MyRecord)) {
        return false;
      }
      return o == this;
    }

    @Override
    public int hashCode() {
      return 0;
    }
  }
  public void test() {
    final SAMFileHeader header = new SAMFileHeader();
    for (int i = 0; i < 100; ++i) {
      header.addSequence(new SAMSequenceRecord("" + i, i + 20));
    }
    final SAMRecord sam1 = new SAMRecord(header);
    sam1.setReadPairedFlag(true);
    sam1.setProperPairFlag(true);
    sam1.setAlignmentStart(23);
    sam1.setMateAlignmentStart(29);
    sam1.setInferredInsertSize(54);
    sam1.setReferenceIndex(42);
    sam1.setMateReferenceIndex(42);
    final MyRecord r = new MyRecord(sam1, 1);

    final SAMRecord sam2 = new SAMRecord(header);
    sam2.setProperPairFlag(true);
    sam2.setReferenceIndex(2);
    sam2.setMateReferenceIndex(12);
    sam2.setMateAlignmentStart(92);
    sam2.setInferredInsertSize(124);
    final MyRecord r2 = new MyRecord(sam2, 12);


    r.setNextInChain(r2);
    assertTrue(r2 == r.chain());

    assertTrue(r.isMated());
    assertFalse(r2.isMated());
    assertEquals(12, r2.getGenome());

    assertEquals(22, r.getStart());
    assertEquals(54, r.getFragmentLength());
    assertEquals(42, r.getMateSequenceId());
    assertEquals(42, r.getSequenceId());

    assertEquals(2, r2.getSequenceId());
    assertEquals(-1, r2.getMateSequenceId());
    assertEquals(-1, r2.getFragmentLength());
  }
}
