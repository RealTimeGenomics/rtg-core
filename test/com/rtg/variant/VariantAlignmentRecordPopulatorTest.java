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

  public void test2() throws IOException {
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
