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
import java.util.List;

import com.rtg.launcher.CommonFlags;
import com.rtg.util.io.InputFileUtils;

/**
 * Class to hold the resulting lists of SAM / calibration files
 */
public class SamCalibrationInputs {
  private final Collection<File> mSamFiles;
  private final Collection<File> mCalibrationFiles;
  /**
   * Constructor for the SAM / calibration file lists
   * @param inputFiles the combined list of input files
   * @param autoload set to true to auto-load calibration files
   * @throws IOException if an IOException occurs
   */
  public SamCalibrationInputs(Collection<File> inputFiles, boolean autoload) throws IOException {
    final List<File> calibrationFiles = new ArrayList<>();
    final List<File> samFiles = new ArrayList<>();
    for (final File f : inputFiles) {
      if (f.getName().endsWith(CommonFlags.RECALIBRATE_EXTENSION)) {
        calibrationFiles.add(f);
      } else {
        samFiles.add(f);
        if (autoload) {
          final File calFile = new File(f.getParent(), f.getName() + CommonFlags.RECALIBRATE_EXTENSION);
          if (calFile.isFile()) {
            calibrationFiles.add(calFile);
          }
        }
      }
    }
    mCalibrationFiles = InputFileUtils.removeRedundantPaths(calibrationFiles);
    mSamFiles = InputFileUtils.removeRedundantPaths(samFiles);
  }

  public Collection<File> getSamFiles() {
    return mSamFiles;
  }

  public Collection<File> getCalibrationFiles() {
    return mCalibrationFiles;
  }
}
