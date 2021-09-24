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
package com.rtg.ngs;

import java.io.File;
import java.util.Arrays;
import java.util.List;

/**
 * Container for record of intermediate files produced during filter concatenation step
 */
public class FilterConcatIntermediateFiles {

  private final List<File> mAlignmentFiles;
  private final List<File> mCalibrationFiles;
  private final List<File> mIndexFiles;

  FilterConcatIntermediateFiles(File[] alignmentFiles, File[] calibrationFiles, File[] indexFiles) {
    mAlignmentFiles = alignmentFiles != null ? Arrays.asList(alignmentFiles) : null;
    mCalibrationFiles = calibrationFiles != null ? Arrays.asList(calibrationFiles) : null;
    mIndexFiles = indexFiles != null ? Arrays.asList(indexFiles) : null;
  }

  public List<File> getAlignmentFiles() {
    return mAlignmentFiles;
  }

  public List<File> getCalibrationFiles() {
    return mCalibrationFiles;
  }

  public List<File> getIndexFiles() {
    return mIndexFiles;
  }
}
