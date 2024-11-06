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

import java.util.Arrays;
import java.util.Collections;

import com.rtg.sam.BadSuperCigarException;
import com.rtg.sam.ReadGroupUtils;
import com.rtg.sam.SamUtils;
import com.rtg.util.StringUtils;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import junit.framework.TestCase;

/**
 */
public class CovariateMachineCycleTest extends TestCase {

  public void testSam() throws BadSuperCigarException {
    final Covariate var = new CovariateMachineCycle(15);
    assertEquals(15, var.size());
    final SAMRecord sam = new SAMRecord(null);

    final Calibrator calwtf = new Calibrator(new Covariate[] {new CovariateReadGroup()}, null);
    final CalibratorCigarParser parser = new CalibratorCigarParser(calwtf) {
      @Override
      public int getReadPosition() {
        return 14;
      }
    };
    assertEquals(14, var.value(sam, parser));
    assertEquals("14", var.valueString(14));
  }

  public void testNameType() {
    final Covariate var = new CovariateMachineCycle(15);
    assertEquals("machinecycle:15", var.name());
    assertEquals(CovariateEnum.MACHINECYCLE, var.getType());
  }

  private byte[] id(final int n) {
    final byte[] q = new byte[n];
    for (int k = 0; k < n; ++k) {
      q[k] = (byte) k;
    }
    return q;
  }

  public void testFR() {
    final Covariate var = new CovariateMachineCycle(0);
    final SAMReadGroupRecord rg = new SAMReadGroupRecord("FR");
    rg.setPlatform("Illumina");
    final SAMFileHeader header = new SAMFileHeader();
    header.addReadGroup(rg);

    final Calibrator cal = new Calibrator(new Covariate[] {var}, null);
    final String refName = "1";
    cal.setTemplate(refName, new byte[10], 10);
    cal.setSequenceLengths(Collections.singletonMap(refName, 10));

    // Left arm forwards
    // |---------->
    //  0123456789
    final SAMRecord sam = new SAMRecord(header);
    sam.setAttribute(ReadGroupUtils.RG_ATTRIBUTE, "FR");
    sam.setReadString("aaaaaaaaaa");
    sam.setBaseQualities(id(10));
    sam.setMappingQuality(40);
    sam.setCigarString("10=");
    sam.setReadPairedFlag(true);
    sam.setFirstOfPairFlag(true);
    sam.setAlignmentStart(1);
    sam.setReferenceName("1");
    cal.processRead(sam);

    // Right arm reverse
    // <----------|
    //  9876543210
    sam.setFirstOfPairFlag(false);
    sam.setSecondOfPairFlag(true);
    sam.setReadNegativeStrandFlag(true);
    sam.setBaseQualities(new byte[] {9, 8, 7, 6, 5, 4, 3, 2, 1, 0});
    cal.processRead(sam);

    int c = 0;
    for (final CalibrationStats stats : cal.mStats) {
      if (stats != null) {
        assertEquals(2, stats.getEqual());
        assertEquals(0, stats.getDifferent());
        assertEquals(0, stats.getDeleted());
        assertEquals(0, stats.getInserted());
        ++c;
      }
    }
    assertEquals(10, c);

//    final ByteArrayOutputStream bos = new ByteArrayOutputStream();
//    cal.writeToStream(bos);
//    System.err.println("" + bos.toString());
  }

  public void testFF() {
    final Covariate var = new CovariateMachineCycle(0);
    final SAMReadGroupRecord rg = new SAMReadGroupRecord("FF");
    rg.setPlatform("Complete");
    final SAMFileHeader header = new SAMFileHeader();
    header.addReadGroup(rg);

    final Calibrator cal = new Calibrator(new Covariate[] {var}, null);
    final String refName = "1";
    cal.setTemplate(refName, new byte[10], 10);
    cal.setSequenceLengths(Collections.singletonMap(refName, 10));

    // Left arm forwards
    // |---------->
    //  0123456789
    final SAMRecord sam = new SAMRecord(header);
    sam.setAttribute(ReadGroupUtils.RG_ATTRIBUTE, "FF");
    sam.setReadString("aaaaaaaaaa");
    sam.setBaseQualities(id(10));
    sam.setMappingQuality(40);
    sam.setCigarString("10=");
    sam.setReadPairedFlag(true);
    sam.setFirstOfPairFlag(true);
    sam.setAlignmentStart(1);
    sam.setReferenceName("1");
    cal.processRead(sam);

    // Right arm forwards
    // |---------->
    //  0123456789
    sam.setFirstOfPairFlag(false);
    sam.setSecondOfPairFlag(true);
    cal.processRead(sam);

    int c = 0;
    for (final CalibrationStats stats : cal.mStats) {
      if (stats != null) {
        assertEquals(2, stats.getEqual());
        assertEquals(0, stats.getDifferent());
        assertEquals(0, stats.getDeleted());
        assertEquals(0, stats.getInserted());
        ++c;
      }
    }
    assertEquals(10, c);
  }

  public void testCGv1() {
    final Covariate var = new CovariateMachineCycle(0);
    final SAMReadGroupRecord rg = new SAMReadGroupRecord("FF");
    rg.setPlatform("COMPLETE");
    final SAMFileHeader header = new SAMFileHeader();
    header.addReadGroup(rg);

    final Calibrator cal = new Calibrator(new Covariate[] {var}, null);
    final byte[] template = new byte[50];
    Arrays.fill(template, (byte) 1);
    final String refName = "1";
    cal.setTemplate(refName, template, 50);
    cal.setSequenceLengths(Collections.singletonMap(refName, 10));

    // Left arm forwards
    // |----------------------------------->
    //  01234567890123456789
    //                 012345678901234
    final SAMRecord sam = new SAMRecord(header);
    sam.setAttribute(ReadGroupUtils.RG_ATTRIBUTE, "FF");
    sam.setReadString(StringUtils.repeat("A", 30));
    sam.setBaseQualities(id(30));
    sam.setMappingQuality(40);
    sam.setAttribute(SamUtils.CG_SUPER_CIGAR, "20=5B15=");
    sam.setAttribute(SamUtils.CG_READ_DELTA, "AAAAA");
    sam.setAttribute(SamUtils.CG_SUPER_CIGAR_OVERLAP_QUALITY, "+++++");
    sam.setReadPairedFlag(true);
    sam.setFirstOfPairFlag(true);
    sam.setAlignmentStart(1);
    sam.setReferenceName("1");
    cal.processRead(sam);

    // Right arm forwards
    // |----------------------------------->
    //  01234567890123456789
    //                 012345678901234
    sam.setFirstOfPairFlag(false);
    sam.setSecondOfPairFlag(true);
    cal.processRead(sam);

    int c = 0;
    for (final CalibrationStats stats : cal.mStats) {
      if (stats != null) {
        assertEquals(2, stats.getEqual());
        assertEquals(0, stats.getDifferent());
        assertEquals(0, stats.getDeleted());
        assertEquals(0, stats.getInserted());
        ++c;
      }
    }
    assertEquals(35, c);
  }

}
