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

package com.rtg.variant;

import java.io.File;

import com.rtg.calibrate.CalibrationStats;
import com.rtg.calibrate.Calibrator;
import com.rtg.calibrate.Calibrator.QuerySpec;
import com.rtg.calibrate.StatsProcessor;
import com.rtg.util.io.FileUtils;
import com.rtg.util.test.FileHelper;

import junit.framework.TestCase;

/**
 */
public class BaseQualityPhredScalerTest extends TestCase {


  public void testABunch() throws Exception {

    final File dir = FileUtils.createTempDir("calMEP", "test");
    try {

      final File cal = FileHelper.resourceToFile("com/rtg/variant/resources/test.calibration", new File(dir, "test.cal"));
      final Calibrator.QuerySpec[] thisQuery = new QuerySpec[1];
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
      assertEquals(14, bqps.getPhred((char) (0 + '!'), 0));
      assertEquals(16, bqps.getPhred((char) (1 + '!'), 0));
      assertEquals(15, bqps.getPhred((char) (2 + '!'), 0));
      assertEquals(15, bqps.getPhred((char) (200 + '!'), 0));
    } finally {
      FileHelper.deleteAll(dir);
    }
  }
}
