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

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;

import com.rtg.launcher.CommonFlags;
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
      final File existingCalFile = new File(tempDir, "samfile0.sam.gz" + CommonFlags.RECALIBRATE_EXTENSION);
      inputFiles.add(existingCalFile);
      for (int i = 0; i < 10; ++i) {
        final File samFile = new File(tempDir, "samfile" + i  + ".sam.gz");
        final File calFile = new File(samFile.getPath() + CommonFlags.RECALIBRATE_EXTENSION);
        final File bamFile = new File(tempDir, "bamfile" + i + ".bam");
        final File bamCalFile = new File(tempDir, "bamcalfile" + i + CommonFlags.RECALIBRATE_EXTENSION);
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
