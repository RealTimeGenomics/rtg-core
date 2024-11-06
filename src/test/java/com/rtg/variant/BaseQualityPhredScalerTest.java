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

package com.rtg.variant;

import java.io.File;

import com.rtg.calibrate.CalibrationStats;
import com.rtg.calibrate.Calibrator;
import com.rtg.calibrate.Calibrator.QuerySpec;
import com.rtg.calibrate.StatsProcessor;
import com.rtg.reader.Arm;
import com.rtg.util.io.TestDirectory;
import com.rtg.util.test.FileHelper;

import junit.framework.TestCase;

/**
 */
public class BaseQualityPhredScalerTest extends TestCase {


  public void testABunch() throws Exception {
    try (final TestDirectory dir = new TestDirectory("calMEP")) {
      final File cal = FileHelper.resourceToFile("com/rtg/variant/resources/test.calibration", new File(dir, "test.cal"));
      final QuerySpec[] thisQuery = new QuerySpec[1];
      final Calibrator c = new Calibrator(Calibrator.getCovariateSet(cal), null) {

        @Override
        public void processStats(StatsProcessor proc, QuerySpec query) {
          assertEquals(thisQuery[0], query);
          final int[] covariateValues = {0, 0, 0, 0};
          final CalibrationStats stats = new CalibrationStats(null);
          stats.seenMnp(3);
          stats.seenEquals(70);
          proc.process(covariateValues, stats);
          final int[] covariateValues2 = {1, 1, 1, 1};
          stats.seenMnp(1);
          stats.seenEquals(92);
          proc.process(covariateValues2, stats);
        }

      };

      thisQuery[0] = c.initQuery();
      final BaseQualityPhredScaler bqps = new BaseQualityPhredScaler(c, thisQuery[0]);
      assertEquals(14, bqps.getScaledPhred((byte) 0, 0, Arm.LEFT));
      assertEquals(16, bqps.getScaledPhred((byte) 1, 0, Arm.LEFT));
      assertEquals(15, bqps.getScaledPhred((byte) 2, 0, Arm.LEFT));
      assertEquals(15, bqps.getScaledPhred((byte) 200, 0, Arm.LEFT));
    }
  }
}
