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
import java.util.HashMap;
import java.util.Map;

import com.rtg.mode.SequenceType;
import com.rtg.reader.MockSequencesReader;
import com.rtg.reference.ReferenceGenome;
import com.rtg.reference.Sex;
import com.rtg.util.TestUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.util.test.FileHelper;
import com.rtg.variant.CalibratedPerSequenceThresholdTest;

import junit.framework.TestCase;

/**
 * Test the corresponding class.
 */
public class ChrStatsTest extends TestCase {

  public void check(final Sex sex, final String... expected) throws IOException {
    final MemoryPrintStream mps = new MemoryPrintStream();
    Diagnostic.setLogStream(mps.printStream());

    final File testDir = FileUtils.createTempDir("junit", "test");
    try {
      FileHelper.resourceToFile("com/rtg/calibrate/resources/reference.txt", new File(testDir, ReferenceGenome.REFERENCE_FILE));

      final HashMap<String, String> readToSample = new HashMap<>();
      readToSample.put("some_reads", "sample");

      final Calibrator c = CalibratedPerSequenceThresholdTest.getCalibrator("com/rtg/calibrate/resources/test.calibration");
      final Map<String, Integer> lengths = new HashMap<>();
      for (int k = 0; k < 9; k++) {
        lengths.put("seq" + k, 10000);
      }
      c.setSequenceLengths(lengths);
      final CalibratedPerSequenceExpectedCoverage calibrator = new CalibratedPerSequenceExpectedCoverage(c, lengths, readToSample, null);
      final MockSequencesReader genomeReader = new MockSequencesReader(SequenceType.DNA, 9) {
        @Override
        public File path() {
          return testDir;
        }
      };
      genomeReader.setLengths(new int[]{10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000});
      try {
        final ChrStats cc = new ChrStats(genomeReader);
        cc.chrStatsCheckAndReport(calibrator, "sample", sex);
      } finally {
        Diagnostic.setLogStream();
      }
      //System.out.println(mps.toString());
      TestUtils.containsAll(mps.toString(), expected);
    } finally {
      assertTrue(FileUtils.deleteFiles(testDir));
    }
  }

  public void testFemaleToMale() throws IOException {
    check(Sex.FEMALE,
      "Average coverage across sequence seq7 for sample sample is 30.32",
      "1 of 8 sequences have unexpected coverage level",
      "FEMALE",
      "MALE");
  }

  public void testEither() throws IOException {
    check(Sex.EITHER,
      "Average coverage across sequence seq7 for sample sample is 30.32",
      "EITHER",
      "MALE");
  }
}
