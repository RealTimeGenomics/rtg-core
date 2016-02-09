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
import java.util.Collections;

import com.rtg.sam.BadSuperCigarException;
import com.rtg.sam.ReadGroupUtils;

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

  public void testResize() throws BadSuperCigarException {
    final Covariate var = new CovariateMachineCycle(6);
    assertEquals(6, var.size());
    final SAMRecord sam = new SAMRecord(null);
    sam.setReadBases(new byte[7]);

    final Calibrator calwtf = new Calibrator(new Covariate[] {new CovariateReadGroup()}, null);
    final CalibratorCigarParser parser = new CalibratorCigarParser(calwtf) {
      private int mPos = 3;
      @Override
      public int getReadPosition() {
        mPos *= 2;
        return mPos;
      }
    };
    assertEquals(6, var.value(sam, parser));
    assertEquals(7, var.newSize());
    assertEquals(6, var.size());
    assertTrue(var.sizeChanged());
    sam.setReadBases(new byte[13]);
    assertEquals(12, var.value(sam, parser));
    assertEquals(13, var.newSize());
    assertEquals(6, var.size());
    assertTrue(var.sizeChanged());
    var.resized();
    assertEquals(13, var.size());
  }

  public void testNameType() {
    final Covariate var = new CovariateMachineCycle(15);
    assertEquals("machinecycle:15", var.name());
    assertEquals(CovariateEnum.MACHINECYCLE, var.getType());
  }

  public void testFR() throws BadSuperCigarException, IOException {
    final Covariate var = new CovariateMachineCycle(0);
    final SAMReadGroupRecord rg = new SAMReadGroupRecord("FR");
    rg.setPlatform("Illumina");
    final SAMFileHeader header = new SAMFileHeader();
    header.addReadGroup(rg);

    final Calibrator cal = new Calibrator(new Covariate[] {var, new CovariateBaseQuality()}, null);
    cal.setTemplate(new byte[10], 10);
    cal.setSequenceLengths(Collections.singletonMap("1", 10));

    // Left arm forwards
    // |---------->
    //  0123456789
    final SAMRecord sam = new SAMRecord(header);
    sam.setAttribute(ReadGroupUtils.RG_ATTRIBUTE, "FR");
    sam.setReadString("aaaaaaaaaa");
    sam.setBaseQualities(new byte[] {0, 1, 2, 3, 4, 5, 6, 7, 8, 9});
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
        c++;
      }
    }
    assertEquals(10, c);

//    final ByteArrayOutputStream bos = new ByteArrayOutputStream();
//    cal.writeToStream(bos);
//    System.err.println("" + bos.toString());
  }

  public void testFF() throws BadSuperCigarException, IOException {
    final Covariate var = new CovariateMachineCycle(0);
    final SAMReadGroupRecord rg = new SAMReadGroupRecord("FF");
    rg.setPlatform("Complete");
    final SAMFileHeader header = new SAMFileHeader();
    header.addReadGroup(rg);

    final Calibrator cal = new Calibrator(new Covariate[] {var, new CovariateBaseQuality()}, null);
    cal.setTemplate(new byte[10], 10);
    cal.setSequenceLengths(Collections.singletonMap("1", 10));

    // Left arm forwards
    // |---------->
    //  0123456789
    final SAMRecord sam = new SAMRecord(header);
    sam.setAttribute(ReadGroupUtils.RG_ATTRIBUTE, "FF");
    sam.setReadString("aaaaaaaaaa");
    sam.setBaseQualities(new byte[] {0, 1, 2, 3, 4, 5, 6, 7, 8, 9});
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
    sam.setReadNegativeStrandFlag(true);
    cal.processRead(sam);

    int c = 0;
    for (final CalibrationStats stats : cal.mStats) {
      if (stats != null) {
        assertEquals(2, stats.getEqual());
        assertEquals(0, stats.getDifferent());
        assertEquals(0, stats.getDeleted());
        assertEquals(0, stats.getInserted());
        c++;
      }
    }
    assertEquals(10, c);
  }

}
