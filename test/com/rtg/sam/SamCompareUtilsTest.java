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

import java.util.Arrays;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.SAMSequenceRecord;

import junit.framework.TestCase;

/**
 */
public class SamCompareUtilsTest extends TestCase {

  public void test() {
    final SAMFileHeader header = new SAMFileHeader();
    header.setSequenceDictionary(new SAMSequenceDictionary(Arrays.asList(new SAMSequenceRecord("raga", 100), new SAMSequenceRecord("yaga", 100), new SAMSequenceRecord("zaga", 100))));
    SAMRecord rec1 = new SAMRecord(header);
    rec1.setReferenceIndex(1);
    SAMRecord rec2 = new SAMRecord(header);
    rec2.setReferenceIndex(2);
    assertEquals(-1, SamCompareUtils.compareSamRecords(rec1, rec2));
    assertEquals(1, SamCompareUtils.compareSamRecords(rec2, rec1));
    rec1.setReferenceIndex(2);
    rec1.setAlignmentStart(50);
    rec2.setAlignmentStart(25);
    assertEquals(1, SamCompareUtils.compareSamRecords(rec1, rec2));
    assertEquals(-1, SamCompareUtils.compareSamRecords(rec2, rec1));
    rec1.setReadPairedFlag(true);
    rec2.setReadPairedFlag(true);
    rec1.setProperPairFlag(true);
    rec2.setProperPairFlag(false);
    rec1.setAlignmentStart(25);
    assertEquals(-1, SamCompareUtils.compareSamRecords(rec1, rec2));
    assertEquals(1, SamCompareUtils.compareSamRecords(rec2, rec1));
    rec2.setProperPairFlag(true);
    rec1.setReadUnmappedFlag(true);
    assertEquals(1, SamCompareUtils.compareSamRecords(rec1, rec2));
    assertEquals(-1, SamCompareUtils.compareSamRecords(rec2, rec1));
    rec2.setReadUnmappedFlag(true);
    assertEquals(0, SamCompareUtils.compareSamRecords(rec1, rec2));
    assertEquals(0, SamCompareUtils.compareSamRecords(rec2, rec1));
  }
}
