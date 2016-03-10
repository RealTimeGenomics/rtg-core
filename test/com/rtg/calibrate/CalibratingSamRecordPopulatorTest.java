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

import java.io.IOException;

import com.rtg.mode.SequenceType;
import com.rtg.reader.MockSequencesReader;
import com.rtg.sam.ReadGroupUtils;

import htsjdk.samtools.SAMRecord;
import junit.framework.TestCase;

/**
 */
public class CalibratingSamRecordPopulatorTest extends TestCase {

  public void test() throws IOException {
    final MockSequencesReader template = new MockSequencesReader(SequenceType.DNA, 5);
    final SAMRecord rec = new SAMRecord(null);
    rec.setAttribute(ReadGroupUtils.RG_ATTRIBUTE, "test");
    rec.setReferenceName(template.name(0));
    final Calibrator cal = new Calibrator(new Covariate[0], null);
    final CalibratingSamRecordPopulator pop = new CalibratingSamRecordPopulator(cal, template, false);
    final SAMRecord r = pop.populate(rec);
    assertEquals(rec, r);
    assertEquals(cal, pop.calibrator());
    rec.setReadUnmappedFlag(true);
    rec.setReferenceName("*");
    final SAMRecord r2 = pop.populate(rec);
    assertEquals(rec, r2);
    assertEquals(cal, pop.calibrator());
  }
}
