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
import java.util.HashMap;
import java.util.Map;

import com.rtg.AbstractTest;
import com.rtg.mode.SequenceType;
import com.rtg.reader.MockSequencesReader;
import com.rtg.reference.ReferenceGenome;
import com.rtg.reference.Sex;
import com.rtg.util.TestUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.util.io.TestDirectory;
import com.rtg.util.test.FileHelper;
import com.rtg.variant.CalibratedPerSequenceThresholdTest;

/**
 * Test the corresponding class.
 */
public class ChrStatsTest extends AbstractTest {

  public void check(final Sex sex, final String... expected) throws IOException {
    try (TestDirectory testDir = new TestDirectory("chrstats")) {
      FileHelper.resourceToFile("com/rtg/calibrate/resources/reference.txt", new File(testDir, ReferenceGenome.REFERENCE_FILE));

      final HashMap<String, String> readToSample = new HashMap<>();
      readToSample.put("some_reads", "sample");

      final Calibrator c = CalibratedPerSequenceThresholdTest.getCalibrator("com/rtg/calibrate/resources/test.calibration");
      final Map<String, Integer> lengths = new HashMap<>();
      for (int k = 0; k < 9; ++k) {
        lengths.put("seq" + k, 10000);
      }
      c.setSequenceLengths(lengths);
      final MockSequencesReader genomeReader = new MockSequencesReader(SequenceType.DNA, 9) {
        @Override
        public File path() {
          return testDir;
        }
      };
      genomeReader.setLengths(new int[]{10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000});

      final MemoryPrintStream mps = new MemoryPrintStream();
      Diagnostic.setLogStream(mps.printStream());
      final CalibratedPerSequenceExpectedCoverage calibrator = new CalibratedPerSequenceExpectedCoverage(c, lengths, readToSample, null);
      final ChrStats cc = new ChrStats(genomeReader);
      cc.chrStatsCheckAndReport(calibrator, "sample", sex);

      //System.out.println(mps.toString());
      TestUtils.containsAll(mps.toString(), expected);
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
