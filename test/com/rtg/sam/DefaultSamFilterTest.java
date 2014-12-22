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
import java.util.ArrayList;
import java.util.Collection;
import java.util.NoSuchElementException;

import com.rtg.util.IntegerOrPercentage;
import com.rtg.util.SingletonPopulatorFactory;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.TestDirectory;
import com.rtg.variant.SamRecordPopulator;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMRecord;

import junit.framework.TestCase;

/**
 */
public class DefaultSamFilterTest extends TestCase {


  public static RecordIterator<SAMRecord> getIterator(File dir, String resource) throws IOException {
    return getIterator(dir, SamFilterParams.builder().create(), resource, false);
  }

  public static RecordIterator<SAMRecord> getIterator(File dir, SamFilterParams filterParams, String resource, boolean indexed) throws IOException {
    final File samFile = new File(dir, resource.substring(resource.lastIndexOf("/") + 1));
    FileUtils.copyResource(resource, samFile);
    if (indexed) {
      final String ext = resource.endsWith(".bam") ? ".bai" : ".tbi";
      FileUtils.copyResource(resource + ext, new File(samFile.getParentFile(), samFile.getName() + ext));
    }
    final Collection<File> coll = new ArrayList<>();
    coll.add(samFile);
    final SamReadingContext c = new SamReadingContext(coll, 1, filterParams, SamUtils.getUberHeader(coll));
    return new ThreadedMultifileIterator<>(c, new SingletonPopulatorFactory<>(new SamRecordPopulator()));
  }

  public void check(SamFilterParams params, int expectedCount, int expectedFiltered, boolean indexed) throws IOException {
    try (TestDirectory dir = new TestDirectory("cnvproducttask")) {
      final RecordIterator<SAMRecord> multiIt = getIterator(dir, params, "com/rtg/variant/cnv/resources/testFilter" + (indexed ? ".bam" : ".sam"), indexed);
      try (SamFilterIterator it = new SamFilterIterator(multiIt, new DefaultSamFilter(params))) {
        int count = 0;
        while (it.hasNext()) {
          it.next();
          count++;
        }
        assertEquals(expectedCount, count);
        assertEquals(expectedFiltered, it.getFilteredRecordsCount());
        assertEquals(false, it.hasNext());
        try {
          it.next();
          fail();
        } catch (final NoSuchElementException e) {
          //ignore
        }
      }
    }
  }

  public void testFilterRecordByFlags() {
    final SamFilterParams.SamFilterParamsBuilder builder = SamFilterParams.builder();
    final SAMRecord rec = new SAMRecord(new SAMFileHeader()); // Not unmapped but alignment position == 0
    rec.setReadUnmappedFlag(true);      // Unmapped with alignment position == 0
    assertTrue(DefaultSamFilter.acceptRecord(builder.create(), rec));
    builder.excludeUnmapped(true);
    assertFalse(DefaultSamFilter.acceptRecord(builder.create(), rec));

    rec.setReadUnmappedFlag(false);      // Mapped with alignment position == 10
    rec.setAlignmentStart(10);
    assertTrue(DefaultSamFilter.acceptRecord(builder.create(), rec));

    rec.setDuplicateReadFlag(true);      // Now a duplicate
    assertTrue(DefaultSamFilter.acceptRecord(builder.create(), rec));
    builder.excludeDuplicates(true);
    assertFalse(DefaultSamFilter.acceptRecord(builder.create(), rec));
    builder.excludeDuplicates(false);

    rec.setReadPairedFlag(true); // Now paired-end

    builder.excludeUnmated(true);
    assertFalse(DefaultSamFilter.acceptRecord(builder.create(), rec));
    rec.setProperPairFlag(true); // Now properly paired (i.e. no longer unmated
    assertTrue(DefaultSamFilter.acceptRecord(builder.create(), rec));

    builder.excludeMated(true);
    assertFalse(DefaultSamFilter.acceptRecord(builder.create(), rec));

  }

  public void testFilterUnplaced() {
    final SamFilterParams.SamFilterParamsBuilder builder = SamFilterParams.builder().excludeUnplaced(true);
    final SAMRecord rec = new SAMRecord(new SAMFileHeader()); // Not unmapped but alignment position == 0
    assertFalse(DefaultSamFilter.acceptRecord(builder.create(), rec));
    rec.setReadUnmappedFlag(true);      // Unmapped with alignment position == 0
    assertFalse(DefaultSamFilter.acceptRecord(builder.create(), rec));
  }

  public void testFilterIterator() throws IOException {
    check(SamFilterParams.builder().create(), 22, 0, false);
    check(SamFilterParams.builder().maxAlignmentCount(-1).maxUnmatedAlignmentScore(new IntegerOrPercentage(0)).create(), 12, 10, false);
    check(SamFilterParams.builder().maxAlignmentCount(0).maxMatedAlignmentScore(new IntegerOrPercentage(10)).create(), 2, 20, false);
    check(SamFilterParams.builder().maxAlignmentCount(1).maxMatedAlignmentScore(null).create(), 22, 0, false);
    check(SamFilterParams.builder().minMapQ(3).excludeUnmapped(true).create(), 17, 5, false);

    // These ones use a restriction which will remove some records prior to being seen by DefaultSamFilter

    try {
      check(SamFilterParams.builder().maxAlignmentCount(-1).maxUnmatedAlignmentScore(new IntegerOrPercentage(0)).restriction("simulatedSequence1").create(), 12, 10, false);
      fail();
    } catch (NoTalkbackSlimException e) {
      assertTrue(e.getMessage().contains("is not indexed"));
    }
    check(SamFilterParams.builder().maxAlignmentCount(-1).maxUnmatedAlignmentScore(new IntegerOrPercentage(0)).restriction("simulatedSequence1").create(), 6, 5, true);
    check(SamFilterParams.builder().maxAlignmentCount(-1).restriction("simulatedSequence1:90-145").create(), 6, 0, true);
  }

}
