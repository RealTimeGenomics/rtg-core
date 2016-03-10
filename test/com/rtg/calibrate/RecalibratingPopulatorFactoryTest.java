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

import com.rtg.mode.SequenceType;
import com.rtg.reader.MockSequencesReader;
import com.rtg.reader.SequencesReader;
import com.rtg.util.Populator;

import htsjdk.samtools.SAMRecord;

import junit.framework.TestCase;

/**
 * Tests the corresponding class.
 *
 */
public class RecalibratingPopulatorFactoryTest extends TestCase {

  public void test() throws IOException {
    final Covariate[] cov = {new CovariateSequence()};
    final SequencesReader t = new MockSequencesReader(SequenceType.DNA, 5);
    final RecalibratingPopulatorFactory f = new RecalibratingPopulatorFactory(cov, null, t);
    final Populator<SAMRecord> p1 = f.populator();
    assertNotSame(p1, f.populator());
    assertNotNull(f.mergedCalibrator());
  }
}
