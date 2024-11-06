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
package com.rtg.variant;

import java.io.File;
import java.io.IOException;

import com.rtg.sam.SamUtils;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.io.TestDirectory;
import com.rtg.util.test.FileHelper;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import junit.framework.TestCase;

/**
 */
public class VariantAlignmentRecordPopulatorTest extends TestCase {

  public void test() {
    final SAMRecord rec = new SAMRecord(new SAMFileHeader());
    rec.setAlignmentStart(42);
    final VariantAlignmentRecordPopulator pop = new VariantAlignmentRecordPopulator(new DefaultMachineErrorChooser(), 0);
    final VariantAlignmentRecord r = pop.populate(rec);
    assertEquals(41, r.getStart());
  }

  public void test2() {
    final VariantAlignmentRecordPopulator pop = new VariantAlignmentRecordPopulator(new DefaultMachineErrorChooser(), 0, "sample1", "sample2");
    final SAMRecord rec = new SAMRecord(new SAMFileHeader());
    rec.setAlignmentStart(42);
    try {
      pop.populate(rec);
      fail("Expected decent exception");
    } catch (NoTalkbackSlimException e) {
      assertTrue(e.getMessage().contains("Encountered a SAM record with no read group information"));
    }
  }

  private static final String SAM = ""
    + "@HD\tVN:1.3\tSO:coordinate\n"
    + "@SQ\tSN:chr13\tLN:285\n"
    + "@RG\tID:ERR000473\tPL:ILLUMINA\n"
    + "4856309\t137\tchr13\t1\t37\t45M\t*\t0\t0\tTCCCAGCTACTCAGGAGTCTGAGGTGGGAGAATTGCTTGATGCCA\tDDB<CCEEDEEHEDGEGEFHCEEF?AFCF??DBA=.=C</:04,7\tRG:Z:ERR000473\n";

  public void test3() throws IOException {
    final VariantAlignmentRecordPopulator pop = new VariantAlignmentRecordPopulator(new DefaultMachineErrorChooser(), 0, "sample1", "sample2");
    try (TestDirectory td = new TestDirectory()) {
      final File samfile = new File(td, "test.sam.gz");
      FileHelper.stringToGzFile(SAM, samfile);
      try (SamReader sr = SamUtils.makeSamReader(samfile)) {
        final SAMRecordIterator iterator = sr.iterator();
        iterator.hasNext();
        final SAMRecord rec = iterator.next();
        try {
          pop.populate(rec);
          fail("Expected decent exception");
        } catch (NoTalkbackSlimException e) {
          assertTrue(e.getMessage(), e.getMessage().contains("Could not determine sample"));
        }
      }
    }
  }
}
