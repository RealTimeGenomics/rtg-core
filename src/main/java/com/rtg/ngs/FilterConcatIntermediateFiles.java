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
