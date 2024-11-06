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

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import com.rtg.reader.SequencesReader;
import com.rtg.util.Populator;
import com.rtg.util.PopulatorFactory;
import com.rtg.util.intervals.ReferenceRegions;

import htsjdk.samtools.SAMRecord;

/**
 * Populator factory where each SAM record populator includes a calibrator.
 */
public class CalibratingPopulatorFactory implements PopulatorFactory<SAMRecord> {

  private final List<CalibratingSamRecordPopulator> mPopulators = new ArrayList<>();

  private final Covariate[] mCovariates;
  private final ReferenceRegions mRegions;
  private final SequencesReader mTemplate;
  private final Map<String, Integer> mLengths;

  /**
   * Construct a populator factory backed by a singleton.
   * @param covariates array of covariates
   * @param regions regions to calibrate
   * @param template the template
   * @throws IOException if there is a problem reading sequence lengths from the reference
   */
  public CalibratingPopulatorFactory(final Covariate[] covariates, final ReferenceRegions regions, final SequencesReader template) throws IOException {
    mCovariates = covariates;
    mRegions = regions;
    mTemplate = template;
    mLengths = (mRegions == null) ? null : Calibrator.getSequenceLengthMap(mTemplate, mRegions);
  }

  @Override
  public Populator<SAMRecord> populator() {
    final Calibrator calibrator = new Calibrator(mCovariates, mRegions);
    if (mLengths != null) {
      calibrator.setSequenceLengths(mLengths);
    }
    final CalibratingSamRecordPopulator pop = new CalibratingSamRecordPopulator(calibrator, mTemplate.copy(), true); // necessary to close these copies?
    mPopulators.add(pop);
    return pop;
  }

  /**
   * Return a calibrator representing the merged content of all the calibrators constructed
   * over the life of this factory.
   * @throws IOException if read or input format errors occur.
   * @return merge calibrator
   */
  public Calibrator mergedCalibrator() throws IOException {
    final Calibrator merged = new Calibrator(mCovariates, null);
    if (mLengths != null) {
      merged.setSequenceLengths(mLengths);
    }
    for (final CalibratingSamRecordPopulator pop : mPopulators) {
      merged.accumulate(pop.calibrator());
    }
    return merged;
  }
}
