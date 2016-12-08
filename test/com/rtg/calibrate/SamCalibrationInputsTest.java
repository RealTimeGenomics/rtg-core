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

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;

import com.rtg.util.io.FileUtils;
import com.rtg.util.test.FileHelper;

import junit.framework.TestCase;

/**
 */
public class SamCalibrationInputsTest extends TestCase {

  public void testSamCalibrationInputs() throws IOException {
    final File tempDir = FileUtils.createTempDir("VarianceInputs", "Test");
    try {
      final Collection<File> inputFiles = new ArrayList<>();
      final StringBuilder sbSamBam = new StringBuilder();
      final StringBuilder sbCalibration = new StringBuilder();
      sbSamBam.append("[");
      sbCalibration.append("[");
      final File existingCalFile = new File(tempDir, "samfile0.sam.gz" + Recalibrate.EXTENSION);
      inputFiles.add(existingCalFile);
      for (int i = 0; i < 10; ++i) {
        final File samFile = new File(tempDir, "samfile" + i  + ".sam.gz");
        final File calFile = new File(samFile.getPath() + Recalibrate.EXTENSION);
        final File bamFile = new File(tempDir, "bamfile" + i + ".bam");
        final File bamCalFile = new File(tempDir, "bamcalfile" + i + Recalibrate.EXTENSION);
        assertTrue(calFile.createNewFile());
        inputFiles.add(samFile);
        inputFiles.add(samFile);
        sbSamBam.append(samFile.getPath());
        sbSamBam.append(", ");
        sbCalibration.append(calFile.getPath());
        sbCalibration.append(", ");
        inputFiles.add(bamFile);
        sbSamBam.append(bamFile.getPath());
        sbSamBam.append(", ");
        inputFiles.add(bamCalFile);
        sbCalibration.append(bamCalFile.getPath());
        sbCalibration.append(", ");
      }
      sbSamBam.setLength(sbSamBam.length() - 2);
      sbSamBam.append("]");
      sbCalibration.setLength(sbCalibration.length() - 2);
      sbCalibration.append("]");
      final SamCalibrationInputs inputs = new SamCalibrationInputs(inputFiles, true);
      assertEquals(20, inputs.getSamFiles().size());
      assertEquals(20, inputs.getCalibrationFiles().size());
      assertEquals(sbSamBam.toString(), inputs.getSamFiles().toString());
      assertEquals(sbCalibration.toString(), inputs.getCalibrationFiles().toString());
    } finally {
      assertTrue(FileHelper.deleteAll(tempDir));
    }
  }
}
