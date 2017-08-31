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
