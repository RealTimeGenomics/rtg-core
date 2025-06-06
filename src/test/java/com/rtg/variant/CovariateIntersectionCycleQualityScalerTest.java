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
import java.io.IOException;
import java.util.function.Consumer;

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
public class CovariateIntersectionCycleQualityScalerTest extends TestCase {

  public void testInterpolationFullRange() throws Exception {
    try (final TestDirectory dir = new TestDirectory("calMEP")) {
      final File cal = FileHelper.resourceToFile("com/rtg/variant/resources/test4.calibration", new File(dir, "test.cal"));
      final QuerySpec[] thisQuery = new QuerySpec[1];
      final Calibrator c = new MyCalibrator(cal, thisQuery,
        proc -> {
          // in readPos 0 from claimed quality 0-63 use interpolated empirical qualities from 10-20
          calibrate(proc, new int[]{0, 0, 0}, 200, 2000);
          calibrate(proc, new int[]{0, 63, 0}, 20, 2000);
        }
      );

      thisQuery[0] = c.initQuery();
      final CovariateIntersectionCycleQualityScaler bqps = new CovariateIntersectionCycleQualityScaler(c, thisQuery[0]);

      // Read position 0
      assertEquals(10, bqps.getScaledPhred((byte) 0, 0, Arm.LEFT));
      assertEquals(15, bqps.getScaledPhred((byte) 31, 0, Arm.LEFT));
      assertEquals(20, bqps.getScaledPhred((byte) 63, 0, Arm.LEFT));

      // We should return the max qual value position when the qual value is too high
      assertEquals(20, bqps.getScaledPhred((byte) 200, 1, Arm.LEFT));

      // We should return the max read position value for read positions that are too high
      assertEquals(10, bqps.getScaledPhred((byte) 0, 1000, Arm.LEFT));
    }
  }
  public void testInterpolation() throws Exception {
    try (final TestDirectory dir = new TestDirectory("calMEP")) {
      final File cal = FileHelper.resourceToFile("com/rtg/variant/resources/test4.calibration", new File(dir, "test.cal"));
      final QuerySpec[] thisQuery = new QuerySpec[1];
      final Calibrator c = new MyCalibrator(cal, thisQuery,
        proc -> {
          // in readPos 10 from claimed quality 20-30 use interpolated empirical qualities from 20-60
          calibrate(proc, new int[]{0, 20, 10}, 100, 10000);
          calibrate(proc, new int[]{0, 30, 10}, 1, 1_000_000);
        }
      );

      thisQuery[0] = c.initQuery();
      final CovariateIntersectionCycleQualityScaler bqps = new CovariateIntersectionCycleQualityScaler(c, thisQuery[0]);

      // Read position 10
      assertEquals(20, bqps.getScaledPhred((byte) 20, 10, Arm.LEFT));
      assertEquals(40, bqps.getScaledPhred((byte) 25, 10, Arm.LEFT));
      assertEquals(60, bqps.getScaledPhred((byte) 30, 10, Arm.LEFT));
    }
  }

  public void testDefaultQualities() throws Exception {
    try (final TestDirectory dir = new TestDirectory("calMEP")) {
      final File cal = FileHelper.resourceToFile("com/rtg/variant/resources/test4.calibration", new File(dir, "test.cal"));
      final QuerySpec[] thisQuery = new QuerySpec[1];
      // No data implies all read positions should range from 0-63
      final Calibrator c = new MyCalibrator(cal, thisQuery,
        proc -> {
        }
      );

      thisQuery[0] = c.initQuery();
      final CovariateIntersectionCycleQualityScaler bqps = new CovariateIntersectionCycleQualityScaler(c, thisQuery[0]);

      for (int i = 0; i < 64; ++i) {
        for (int j = 0; j < 35; ++j) {
          assertEquals(String.format("readPos:%d, claimed quality: %d", j, i), i, bqps.getScaledPhred((byte) i, j, Arm.LEFT));
        }
      }
    }
  }

  public void testAllPositionQualities() throws Exception {
    try (final TestDirectory dir = new TestDirectory("calMEP")) {
      final File cal = FileHelper.resourceToFile("com/rtg/variant/resources/test4.calibration", new File(dir, "test.cal"));
      final QuerySpec[] thisQuery = new QuerySpec[1];
      // No data implies all read positions should range from 0-63
      final Calibrator c = new MyCalibrator(cal, thisQuery,
        proc -> {
          // None of these are in the same read position so should be no interpolation between per read position quality

          // at qual 20 twice the error rate
          calibrate(proc, new int[]{0, 20, 5}, 20, 1000);
          calibrate(proc, new int[]{0, 20, 10}, 20, 1000);

          // at qual 40 half the error rate
          calibrate(proc, new int[]{0, 40, 15}, 0, 10_000);
          calibrate(proc, new int[]{0, 40, 20}, 1, 10_000);
        }

      );

      thisQuery[0] = c.initQuery();
      final CovariateIntersectionCycleQualityScaler bqps = new CovariateIntersectionCycleQualityScaler(c, thisQuery[0]);

      // at read position 0
      // we should be interpolated from claimed (0-20) observed (0-17)
      assertEquals(0, bqps.getScaledPhred((byte) 0, 0, Arm.LEFT));
      assertEquals(9, bqps.getScaledPhred((byte) 10, 0, Arm.LEFT));
      assertEquals(17, bqps.getScaledPhred((byte) 20, 0, Arm.LEFT));

      //from claimed (20-40) we should observe (17-43)
      assertEquals(30, bqps.getScaledPhred((byte) 30, 0, Arm.LEFT));
      assertEquals(43, bqps.getScaledPhred((byte) 40, 0, Arm.LEFT));
    }
  }

  public void testFromCovariates() throws Exception {
    try (final TestDirectory dir = new TestDirectory("calMEP")) {
      final File cal = FileHelper.resourceToFile("com/rtg/variant/resources/test4.calibration", new File(dir, "test.cal"));
      final QuerySpec[] thisQuery = new QuerySpec[1];
      // No data implies all read positions should range from 0-63
      final Calibrator c = new MyCalibrator(cal, thisQuery,
        proc -> {
          // We place points around a central one. qualities known at * and imputed at +, X is there to shift the global quality
          //   X *
          //
          //   * + *
          //
          //     *


          // at cycle 15 the average quality is 20 (measured at claimed qualities 21, 40
          calibrate(proc, new int[]{0, 10, 15}, 5000, 1_000_000);
          calibrate(proc, new int[]{0, 30, 15}, 15000, 1_000_000);

          // at quality 20 the average quality is 40
          calibrate(proc, new int[]{0, 20, 10}, 50, 1_000_000);
          calibrate(proc, new int[]{0, 20, 20}, 150, 1_000_000);

          // Somewhere else our global quality is worse...
          calibrate(proc, new int[]{0, 10, 20}, 15000, 1_000_000);

        }

      );

      thisQuery[0] = c.initQuery();
      final CovariateIntersectionCycleQualityScaler bqps = new CovariateIntersectionCycleQualityScaler(c, thisQuery[0]);

      assertEquals(0, bqps.getScaledPhred((byte) 0, 0, Arm.LEFT));

      // at read position 0
      // we should be interpolated from claimed (0-20) observed (0-17)
      assertEquals(30, bqps.getScaledPhred((byte) 20, 15, Arm.LEFT));
    }
  }

  private void calibrate(StatsProcessor proc, int[] covariateValues, int mnps, int equals) {
    final CalibrationStats stats = new CalibrationStats(null);
    stats.seenMnp(mnps);
    stats.seenEquals(equals);
    proc.process(covariateValues, stats);
  }

  private class MyCalibrator extends Calibrator {

    private final QuerySpec[] mThisQuery;
    private final Consumer<StatsProcessor> mConsumer;

    MyCalibrator(File cal, QuerySpec[] thisQuery, Consumer<StatsProcessor> consumer) throws IOException {
      super(getCovariateSet(cal), null);
      mThisQuery = thisQuery;
      mConsumer = consumer;
    }

    @Override
    public void processStats(StatsProcessor proc, QuerySpec query) {
      assertEquals(mThisQuery[0], query);
      mConsumer.accept(proc);
    }
  }

}
