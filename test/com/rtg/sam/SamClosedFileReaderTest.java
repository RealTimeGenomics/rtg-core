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

import java.io.File;
import java.io.IOException;
import java.util.Iterator;

import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.intervals.ReferenceRanges;
import com.rtg.util.io.FileUtils;
import com.rtg.util.test.FileHelper;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;

import junit.framework.TestCase;

/**
 * Test class
 */
public class SamClosedFileReaderTest extends TestCase {

  @Override
  protected void setUp() {
    Diagnostic.setLogStream();
  }

  public void testSam() throws IOException {
    final File dir = FileUtils.createTempDir("closedfilereader", "test");
    try {
      final File sam = FileHelper.resourceToFile("com/rtg/variant/bayes/multisample/cancer/resources/test1_normal.sam", new File(dir, "sam.sam"));
      check(sam);
    } finally {
      assertTrue(FileHelper.deleteAll(dir));
    }
  }

  private void check(final File bam) throws IOException {
    final SAMFileReader normalReader = new SAMFileReader(bam);
    final Iterator<SAMRecord> norm = normalReader.iterator();
    try {
      try (final RecordIterator<SAMRecord> closed = new SamClosedFileReader(bam, null, normalReader.getFileHeader())) {
        while (norm.hasNext() && closed.hasNext()) {
          final SAMRecord a = norm.next();
          final SAMRecord b = closed.next();
          assertEquals(a.getSAMString().trim(), b.getSAMString().trim());
        }
        assertEquals(norm.hasNext(), closed.hasNext());
      }
    } finally {
      normalReader.close();
    }
  }

  public void testBam() throws IOException {
    final File dir = FileUtils.createTempDir("closedfilereader", "test");
    try {
      final File sam = FileHelper.resourceToFile("com/rtg/variant/bayes/multisample/cancer/resources/test1_normal.sam", new File(dir, "sam.sam"));
      final File bam = new File(dir, "bam.bam");
      Sam2Bam.convertSamToBam(bam, BamIndexer.indexFileName(bam), sam);
      check(bam);
    } finally {
      assertTrue(FileHelper.deleteAll(dir));
    }
  }

  private void check2(final File sam, File exp) throws IOException {
    try (SAMFileReader normalReader = new SAMFileReader(exp)) {
      final SamRegionRestriction rr = new SamRegionRestriction("simulatedSequence2", 10000, 20000);
      final Iterator<SAMRecord> norm = normalReader.query(rr.getSequenceName(), rr.getStart() + 1, rr.getEnd(), false);
      ReferenceRanges ranges = SamRangeUtils.createExplicitReferenceRange(normalReader.getFileHeader(), rr);
      try (final RecordIterator<SAMRecord> closed = new SamClosedFileReader(sam, ranges, normalReader.getFileHeader())) {
        while (norm.hasNext() && closed.hasNext()) {
          final SAMRecord a = norm.next();
          final SAMRecord b = closed.next();
          assertEquals(a.getSAMString(), b.getSAMString());
        }
        assertEquals(norm.hasNext(), closed.hasNext());
      }
    }
  }

  public void testSamRestrict() throws IOException {
    final File dir = FileUtils.createTempDir("closedfilereader", "test");
    try {
      final File sam = FileHelper.resourceToFile("com/rtg/sam/resources/readerWindow2.sam.gz", new File(dir, "sam.sam.gz"));
      FileHelper.resourceToFile("com/rtg/sam/resources/readerWindow2.sam.gz.tbi", new File(dir, "sam.sam.gz.tbi"));
      final File bam = new File(dir, "bam.bam");
      Sam2Bam.convertSamToBam(bam, BamIndexer.indexFileName(bam), sam);
      check2(sam, bam);
    } finally {
      assertTrue(FileHelper.deleteAll(dir));
    }
  }

  public void testBamRestrict() throws IOException {
    final File dir = FileUtils.createTempDir("closedfilereader", "test");
    try {
      final File sam = FileHelper.resourceToFile("com/rtg/sam/resources/readerWindow2.sam.gz", new File(dir, "sam.sam.gz"));
      final File bam = new File(dir, "bam.bam");
      Sam2Bam.convertSamToBam(bam, BamIndexer.indexFileName(bam), sam);
      check2(bam, bam);
    } finally {
      assertTrue(FileHelper.deleteAll(dir));
    }
  }
}
