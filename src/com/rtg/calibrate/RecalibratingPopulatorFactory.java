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
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import com.rtg.reader.SequencesReader;
import com.rtg.util.Populator;
import com.rtg.util.PopulatorFactory;
import com.rtg.util.intervals.ReferenceRegions;
import com.rtg.variant.RecalibratingSamRecordPopulator;

import htsjdk.samtools.SAMRecord;

/**
 * Populator factory where each SAM record populator includes a calibrator.
 */
public class RecalibratingPopulatorFactory implements PopulatorFactory<SAMRecord> {

  private final List<RecalibratingSamRecordPopulator> mPopulators = new ArrayList<>();

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
  public RecalibratingPopulatorFactory(final Covariate[] covariates, final ReferenceRegions regions, final SequencesReader template) throws IOException {
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
    final RecalibratingSamRecordPopulator pop = new RecalibratingSamRecordPopulator(calibrator, mTemplate.copy(), true); // necessary to close these copies?
    mPopulators.add(pop);
    return pop;
  }

  /**
   * Return a calibrator representing the merged content of all the calibrators constructed
   * over the life of this factory.
   * @return merge calibrator
   */
  public Calibrator mergedCalibrator() {
    final Calibrator merged = new Calibrator(mCovariates, null);
    if (mLengths != null) {
      merged.setSequenceLengths(mLengths);
    }
    for (final RecalibratingSamRecordPopulator pop : mPopulators) {
      merged.accumulate(pop.calibrator());
    }
    return merged;
  }
}
