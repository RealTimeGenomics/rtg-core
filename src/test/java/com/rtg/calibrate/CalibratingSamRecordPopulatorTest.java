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
package com.rtg.calibrate;

import com.rtg.mode.SequenceType;
import com.rtg.reader.MockSequencesReader;
import com.rtg.sam.ReadGroupUtils;

import htsjdk.samtools.SAMRecord;
import junit.framework.TestCase;

/**
 */
public class CalibratingSamRecordPopulatorTest extends TestCase {

  public void test() {
    final MockSequencesReader template = new MockSequencesReader(SequenceType.DNA, 5);
    final SAMRecord rec = new SAMRecord(null);
    rec.setAttribute(ReadGroupUtils.RG_ATTRIBUTE, "test");
    rec.setReferenceName(template.name(0));
    rec.setCigarString("2=");
    rec.setReadString("ac");
    rec.setBaseQualityString("D!");
    rec.setAlignmentStart(1);
    rec.setMappingQuality(1);
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
