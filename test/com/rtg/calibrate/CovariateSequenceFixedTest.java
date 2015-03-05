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
package com.rtg.calibrate;

import java.util.ArrayList;

import com.rtg.sam.SamUtils;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;

import junit.framework.TestCase;

/**
 */
public class CovariateSequenceFixedTest extends TestCase {

  private SAMRecord recSam(String sequence, SAMFileHeader header) {
    final SAMRecord sam = new SAMRecord(header);
    sam.setReferenceName(sequence);
    sam.setCigarString("35=");
    sam.setReadString("acgtacgtggacgtacgtggacgtacgtggttttt");
    sam.setBaseQualityString("012345678901234567890123456789ABCDE");
    sam.setAlignmentStart(1);
    sam.setMappingQuality(1);
    return sam;
  }

  private SAMFileHeader makeHeader() {
    final SAMFileHeader header = new SAMFileHeader();
    final SAMSequenceDictionary dict = header.getSequenceDictionary();
    dict.addSequence(new SAMSequenceRecord("sequence1", 100));
    dict.addSequence(new SAMSequenceRecord("sequence2", 100));
    dict.addSequence(new SAMSequenceRecord("sequence3", 100));
    return header;
  }

  public void testValue() {
    final SAMFileHeader header = makeHeader();
    final CovariateSequenceFixed cs = new CovariateSequenceFixed(SamUtils.getSequenceNames(header));
    final Calibrator cal = new Calibrator(new Covariate[] {cs, new CovariateReadGroup()}, null);
    final CalibratorCigarParser parser = new CalibratorCigarParser(cal);
    assertEquals(3, cs.newSize());
    assertEquals(3, cs.size());

    assertEquals(0, cs.value(recSam("sequence1", header), parser));
    assertEquals(1, cs.value(recSam("sequence2", header), parser));
    assertEquals(2, cs.value(recSam("sequence3", header), parser));
    assertEquals(0, cs.value(recSam("sequence1", header), parser));
    assertEquals(2, cs.value(recSam("sequence3", header), parser));

  }

  public void testValueString() {
    final SAMFileHeader header = makeHeader();
    final CovariateSequenceFixed cs = new CovariateSequenceFixed(SamUtils.getSequenceNames(header));
    final Calibrator cal = new Calibrator(new Covariate[] {cs, new CovariateReadGroup()}, null);
    final CalibratorCigarParser parser = new CalibratorCigarParser(cal);

    cs.value(recSam("sequence1", header), parser);
    cs.value(recSam("sequence2", header), parser);
    cs.value(recSam("sequence3", header), parser);
    assertEquals("sequence1", cs.valueString(0));
    assertEquals("sequence2", cs.valueString(1));
    assertEquals("sequence3", cs.valueString(2));

  }

  public void testParse() {
    final SAMFileHeader header = makeHeader();
    final CovariateSequenceFixed cs = new CovariateSequenceFixed(SamUtils.getSequenceNames(header));
    assertEquals(0, cs.parse("sequence1"));
    try {
      cs.parse("no-such-sequence1");
      fail();
    } catch (final UnsupportedOperationException e) {
      // ok
    }
  }

  public void testGetType() {
    final CovariateSequenceFixed cs = new CovariateSequenceFixed(new ArrayList<String>());
    assertEquals(CovariateEnum.SEQUENCE, cs.getType());
  }
}
