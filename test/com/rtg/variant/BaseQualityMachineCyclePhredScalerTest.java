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
import java.io.IOException;
import java.util.function.Consumer;

import com.rtg.calibrate.CalibrationStats;
import com.rtg.calibrate.Calibrator;
import com.rtg.calibrate.Calibrator.QuerySpec;
import com.rtg.calibrate.StatsProcessor;
import com.rtg.ngs.Arm;
import com.rtg.util.io.FileUtils;
import com.rtg.util.test.FileHelper;

import junit.framework.TestCase;

/**
 */
public class BaseQualityMachineCyclePhredScalerTest extends TestCase {


  public void testABunch() throws Exception {

    final File dir = FileUtils.createTempDir("calMEP", "test");
    try {

      final File cal = FileHelper.resourceToFile("com/rtg/variant/resources/test4.calibration", new File(dir, "test.cal"));
      final Calibrator.QuerySpec[] thisQuery = new QuerySpec[1];
      final Calibrator c = new MyCalibrator(cal, thisQuery,
        proc -> {
          calibrate(proc, new int[]{0, 0, 0, 0}, 300, 7000);
          calibrate(proc, new int[]{1, 1, 1, 1}, 100, 9200);
        }
      );

      thisQuery[0] = c.initQuery();
      final BaseQualityMachineCyclePhredScaler bqps = new BaseQualityMachineCyclePhredScaler(c, thisQuery[0]);
      assertEquals(14, bqps.getScaledPhred((byte) 0, 0, Arm.LEFT));
      assertEquals(20, bqps.getScaledPhred((byte) 1, 0, Arm.LEFT));
      assertEquals(21, bqps.getScaledPhred((byte) 2, 0, Arm.LEFT));
      assertEquals(14, bqps.getScaledPhred((byte) 0, 1, Arm.LEFT));
      assertEquals(20, bqps.getScaledPhred((byte) 1, 1, Arm.LEFT));
      assertEquals(21, bqps.getScaledPhred((byte) 2, 1, Arm.LEFT));

      // We should return the max qual value position when the qual value is too high
      assertEquals(63, bqps.getScaledPhred((byte) 200, 1, Arm.LEFT));

      // We should return the max read position value for read positions that are too high
      assertEquals(14, bqps.getScaledPhred((byte) 0, 1000, Arm.LEFT));
    } finally {
      FileHelper.deleteAll(dir);
    }
  }

  public void testInterpolationFullRange() throws Exception {

    final File dir = FileUtils.createTempDir("calMEP", "test");
    try {

      final File cal = FileHelper.resourceToFile("com/rtg/variant/resources/test4.calibration", new File(dir, "test.cal"));
      final Calibrator.QuerySpec[] thisQuery = new QuerySpec[1];
      final Calibrator c = new MyCalibrator(cal, thisQuery,
        proc -> {
          // in readPos 0 from claimed quality 0-63 use interpolated empirical qualities from 10-20
          calibrate(proc, new int[]{0, 0, 0}, 200, 2000);
          calibrate(proc, new int[]{0, 63, 0}, 20, 2000);
        }
      );

      thisQuery[0] = c.initQuery();
      final BaseQualityMachineCyclePhredScaler bqps = new BaseQualityMachineCyclePhredScaler(c, thisQuery[0]);

      // Read position 0
      assertEquals(10, bqps.getScaledPhred((byte) 0, 0, Arm.LEFT));
      assertEquals(15, bqps.getScaledPhred((byte) 31, 0, Arm.LEFT));
      assertEquals(20, bqps.getScaledPhred((byte) 63, 0, Arm.LEFT));

      // We should return the max qual value position when the qual value is too high
      assertEquals(20, bqps.getScaledPhred((byte) 200, 1, Arm.LEFT));

      // We should return the max read position value for read positions that are too high
      assertEquals(10, bqps.getScaledPhred((byte) 0, 1000, Arm.LEFT));
    } finally {
      FileHelper.deleteAll(dir);
    }
  }
  public void testInterpolation() throws Exception {

    final File dir = FileUtils.createTempDir("calMEP", "test");
    try {

      final File cal = FileHelper.resourceToFile("com/rtg/variant/resources/test4.calibration", new File(dir, "test.cal"));
      final Calibrator.QuerySpec[] thisQuery = new QuerySpec[1];
      final Calibrator c = new MyCalibrator(cal, thisQuery,
        proc -> {
          // in readPos 10 from claimed quality 20-30 use interpolated empirical qualities from 20-60
          calibrate(proc, new int[]{0, 20, 10}, 100, 10_000);
          calibrate(proc, new int[]{0, 30, 10}, 100, 100_000_000);
        }
      );

      thisQuery[0] = c.initQuery();
      final BaseQualityMachineCyclePhredScaler bqps = new BaseQualityMachineCyclePhredScaler(c, thisQuery[0]);

      // Read position 10
      assertEquals(20, bqps.getScaledPhred((byte) 20, 10, Arm.LEFT));
      assertEquals(40, bqps.getScaledPhred((byte) 25, 10, Arm.LEFT));
      assertEquals(60, bqps.getScaledPhred((byte) 30, 10, Arm.LEFT));

    } finally {
      FileHelper.deleteAll(dir);
    }
  }

  public void testDefaultQualities() throws Exception {

    final File dir = FileUtils.createTempDir("calMEP", "test");
    try {

      final File cal = FileHelper.resourceToFile("com/rtg/variant/resources/test4.calibration", new File(dir, "test.cal"));
      final Calibrator.QuerySpec[] thisQuery = new QuerySpec[1];
      // No data implies all read positions should range from 0-63
      final Calibrator c = new MyCalibrator(cal, thisQuery,
        proc -> { }
      );

      thisQuery[0] = c.initQuery();
      final BaseQualityMachineCyclePhredScaler bqps = new BaseQualityMachineCyclePhredScaler(c, thisQuery[0]);

      for (int i = 0; i < 64; i++) {
        for (int j = 0; j < 35; j++) {
          assertEquals(String.format("readPos:%d, claimed quality: %d", j, i), i, bqps.getScaledPhred((byte) i , j, Arm.LEFT));
        }
      }

    } finally {
      FileHelper.deleteAll(dir);
    }
  }

  public void testAllPositionQualities() throws Exception {

    final File dir = FileUtils.createTempDir("calMEP", "test");
    try {

      final File cal = FileHelper.resourceToFile("com/rtg/variant/resources/test4.calibration", new File(dir, "test.cal"));
      final Calibrator.QuerySpec[] thisQuery = new QuerySpec[1];
      // No data implies all read positions should range from 0-63
      final Calibrator c = new MyCalibrator(cal, thisQuery,
        proc -> {
          // None of these are in the same read position so should be no interpolation between per read position quality

          // at qual 20 twice the error rate
          calibrate(proc, new int[]{0, 20, 5}, 20, 1000);
          calibrate(proc, new int[]{0, 20, 10}, 20, 1000);

          // at qual 40 half the error rate
          calibrate(proc, new int[]{0, 40, 15}, 10, 100_000);
          calibrate(proc, new int[]{0, 40, 20}, 10, 100_000);
        }

      );

      thisQuery[0] = c.initQuery();
      final BaseQualityMachineCyclePhredScaler bqps = new BaseQualityMachineCyclePhredScaler(c, thisQuery[0]);

      // at read position 0
      // we should be interpolated from claimed (0-20) observed (0-17)
      assertEquals(0, bqps.getScaledPhred((byte) 0, 0, Arm.LEFT));
      assertEquals(9, bqps.getScaledPhred((byte) 10, 0, Arm.LEFT));
      assertEquals(17, bqps.getScaledPhred((byte) 20, 0, Arm.LEFT));

      //from claimed (20-40) we should observe (17-40)
      assertEquals(29, bqps.getScaledPhred((byte) 30, 0, Arm.LEFT));
      assertEquals(40, bqps.getScaledPhred((byte) 40, 0, Arm.LEFT));

    } finally {
      FileHelper.deleteAll(dir);
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
      super(Calibrator.getCovariateSet(cal), null);
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
