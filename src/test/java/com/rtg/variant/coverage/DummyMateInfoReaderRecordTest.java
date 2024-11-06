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
