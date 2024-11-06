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

import com.rtg.calibrate.Calibrator;
import com.rtg.calibrate.Covariate;
import com.rtg.reader.Arm;
import com.rtg.sam.ReadGroupUtils;
import com.rtg.util.InvalidParamsException;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.io.FileUtils;
import com.rtg.util.test.FileHelper;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import junit.framework.TestCase;

/**
 * Test class
 */
public class CalibratedMachineErrorChooserTest extends TestCase {

  public void testNormal() throws IOException, InvalidParamsException {
    checkCalibration("com/rtg/variant/resources/test2.calibration", 3, 3);
  }

  public void testLong() throws IOException, InvalidParamsException {
    checkCalibration("com/rtg/variant/resources/test3.calibration", 30, 30);
  }

  public void testLongReadPos() throws IOException, InvalidParamsException {
    checkCalibration("com/rtg/variant/resources/test4.calibration", 40, 30);
  }


  Calibrator getCalibrator(File... calibrationFiles) throws IOException {
    final Covariate[] covariates;
    if (calibrationFiles.length > 0) {
      covariates = Calibrator.getCovariateSet(calibrationFiles[0]);
    } else {
      throw new IllegalArgumentException("need calibration files");
    }
    final Calibrator c = new Calibrator(covariates, null);
    for (File calibrationFile : calibrationFiles) {
      c.accumulate(calibrationFile);
    }
    return c;
  }

  public void checkCalibration(String file, int expectedQual0, int expectedQual1) throws IOException, InvalidParamsException {
    Diagnostic.setLogStream();
    final File dir = FileUtils.createTempDir("calMEC", "test");
    try {
      final File cal = FileHelper.resourceToFile(file, new File(dir, "test.cal"));
      final CalibratedMachineErrorChooser c = new CalibratedMachineErrorChooser(getCalibrator(cal));
      final SAMFileHeader fheader = new SAMFileHeader();
      final SAMReadGroupRecord srgr1 = new SAMReadGroupRecord("readgroup1");
      srgr1.setPlatform("ILLUMINA");
      final SAMReadGroupRecord srgr2 = new SAMReadGroupRecord("readgroup2");
      srgr2.setPlatform("ILLUMINA");
      fheader.addReadGroup(srgr1);
      fheader.addReadGroup(srgr2);
      final SAMRecord fake = new SAMRecord(fheader);
      VariantAlignmentRecord var = new VariantAlignmentRecord(fake);
      try {
        c.machineErrors(var.getReadGroup(), var.isReadPaired());
        fail();
      } catch (final NoTalkbackSlimException e) {
        assertTrue(e.getMessage().contains("Read group"));
      }
      fake.setAttribute(ReadGroupUtils.RG_ATTRIBUTE, "readgroup1");  // known
      var = new VariantAlignmentRecord(fake);
      assertEquals(expectedQual0, c.machineErrors(var.getReadGroup(), var.isReadPaired()).getScaledPhredFromAscii((char) ('!' + 2), 0, Arm.LEFT));
      assertEquals(expectedQual1, c.machineErrors(var.getReadGroup(), var.isReadPaired()).getScaledPhredFromAscii((char) ('!' + 2), 1, Arm.LEFT));
      fake.setAttribute(ReadGroupUtils.RG_ATTRIBUTE, "readgroup2");  // unknown
      try {
        var = new VariantAlignmentRecord(fake);
        c.machineErrors(var.getReadGroup(), var.isReadPaired());
        fail();
      } catch (final NoTalkbackSlimException e) {
        assertEquals("No calibration data supplied for read group: readgroup2", e.getMessage());
      }
    } finally {
      assertTrue(FileHelper.deleteAll(dir));
    }
  }
}
