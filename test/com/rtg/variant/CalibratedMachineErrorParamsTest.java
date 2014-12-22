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
import java.util.HashMap;
import java.util.Map;

import com.rtg.calibrate.CalibrationStats;
import com.rtg.calibrate.Calibrator;
import com.rtg.calibrate.Covariate;
import com.rtg.calibrate.CovariateBaseQuality;
import com.rtg.calibrate.CovariateEnum;
import com.rtg.util.Histogram;
import com.rtg.util.io.FileUtils;
import com.rtg.util.machine.MachineType;
import com.rtg.util.test.FileHelper;

import junit.framework.TestCase;

/**
 * Test class
 */
public class CalibratedMachineErrorParamsTest extends TestCase {

  public void testDistributions() {
    final Histogram gap = new Histogram();
    final Histogram lap = new Histogram();

    //0 is 7, but 0 position is never stored in file
    gap.increment(1, 6);
    gap.increment(2, 5);
    gap.increment(3, 4);
    gap.increment(4, 4);
    gap.increment(5, 5);
    gap.increment(6, 6);
    gap.increment(7, 7);

    //0 is 4
    lap.increment(1, 7);
    lap.increment(2, 6);
    lap.increment(3, 5);

    double[][] dists = CalibratedMachineErrorParams.histogramsToCgDists(gap, lap);
    //gap dists are length 4 + index
    assertEquals((double) 4 / 22, dists[0][0]);
    assertEquals((double) 5 / 22, dists[0][1]);
    assertEquals((double) 6 / 22, dists[0][2]);
    assertEquals((double) 7 / 22, dists[0][3]);
    assertEquals((double) 0 / 22, dists[0][4]);
    //small gaps are length == index
    assertEquals((double) 7 / 22, dists[1][0]);
    assertEquals((double) 6 / 22, dists[1][1]);
    assertEquals((double) 5 / 22, dists[1][2]);
    assertEquals((double) 4 / 22, dists[1][3]);
    //overlaps are back to front (position 0 is for length 4, position 4 is for length 0);
    assertEquals((double) 0 / 22, dists[2][0]);
    assertEquals((double) 5 / 22, dists[2][1]);
    assertEquals((double) 6 / 22, dists[2][2]);
    assertEquals((double) 7 / 22, dists[2][3]);
    assertEquals((double) 4 / 22, dists[2][4]);
  }

  public void testRates() {
    final Calibrator c = new Calibrator(new Covariate[] {new CovariateBaseQuality()}, null) {
      @Override
      public CalibrationStats getSums(CovariateEnum covariate, String covariateValue) {
        final CalibrationStats stats = new CalibrationStats(null);
        stats.seenMnp(40);
        stats.seenInsert(30);
        stats.seenDelete(20);
        stats.seenEquals(100);
        return stats;
      }
      @Override
      public Histogram getHistogram(String type, String readGroup) {
        final Histogram hist = new Histogram();
        hist.increment(1);
        return hist;
      }
    };
    final CalibratedMachineErrorParams cp = new CalibratedMachineErrorParams(MachineType.ILLUMINA_PE, c, "unknown");
    assertEquals((double) 30 / 141, cp.errorInsBaseRate());
    assertEquals((double) 20 / 141, cp.errorDelBaseRate());
    assertEquals((double) 40 / 141, cp.errorSnpRate());
    //distribution set such that event -> length = 1 when no distribution available
    assertEquals((double) 30 / 141, cp.errorInsEventRate());
    assertEquals((double) 20 / 141, cp.errorDelEventRate());
    assertEquals((double) 40 / 141, cp.errorMnpEventRate());
  }
  static final Map<Integer, Integer> TEST_CALIBRATION_QUALITIES = new HashMap<>();
  static {
    final int totalEqPlusMm = 214733;
    final int totalMm = 4673;
    double errorRate = (double) totalMm / (double) totalEqPlusMm;
    TEST_CALIBRATION_QUALITIES.put(2, CalibratedMachineErrorParams.countsToEmpiricalQuality(270, 270 + 467, errorRate));
    TEST_CALIBRATION_QUALITIES.put(20, CalibratedMachineErrorParams.countsToEmpiricalQuality(289, 289 + 3434, errorRate));
    TEST_CALIBRATION_QUALITIES.put(30, CalibratedMachineErrorParams.countsToEmpiricalQuality(86, 86 + 16424, errorRate));
    TEST_CALIBRATION_QUALITIES.put(32, CalibratedMachineErrorParams.countsToEmpiricalQuality(248, 248 + 52395, errorRate));

    TEST_CALIBRATION_QUALITIES.put(36, CalibratedMachineErrorParams.countsToEmpiricalQuality(1, 1 + 321, errorRate));
    TEST_CALIBRATION_QUALITIES.put(63, CalibratedMachineErrorParams.countsToEmpiricalQuality(0, 0, errorRate));

  }

  public void testMinimumPhred() {
    assertEquals(2, CalibratedMachineErrorParams.countsToEmpiricalQuality(5000, 6000, 0.99));
    assertEquals(2, CalibratedMachineErrorParams.countsToEmpiricalQuality(5000, 6000, 0.001));
    assertEquals(2, CalibratedMachineErrorParams.countsToEmpiricalQuality(5000, 6000, 0.001));
    assertEquals(2, CalibratedMachineErrorParams.countsToEmpiricalQuality(5630, 10000, 0.001));
    assertEquals(3, CalibratedMachineErrorParams.countsToEmpiricalQuality(5620, 10000, 0.001));
    assertEquals(10, CalibratedMachineErrorParams.countsToEmpiricalQuality(1000, 10000, 0.001));
  }

  public void testCurve() throws IOException {
    final File dir = FileUtils.createTempDir("calMEP", "test");
    try {
      final File cal = FileHelper.resourceToFile("com/rtg/variant/resources/test.calibration", new File(dir, "test.cal"));
      final Calibrator c = new Calibrator(Calibrator.getCovariateSet(cal), null);
      c.accumulate(cal);
      final CalibratedMachineErrorParams cp = new CalibratedMachineErrorParams(MachineType.ILLUMINA_PE, c, "readgroup1");

      assertEquals((int) TEST_CALIBRATION_QUALITIES.get(2), cp.getPhred((char) ('!' + 2), 1));
      assertEquals((int) TEST_CALIBRATION_QUALITIES.get(20), cp.getPhred((char) ('!' + 20), 1));
      assertEquals((int) TEST_CALIBRATION_QUALITIES.get(30), cp.getPhred((char) ('!' + 30), 1));
      assertEquals((int) TEST_CALIBRATION_QUALITIES.get(32), cp.getPhred((char) ('!' + 32), 1));
      assertEquals((int) TEST_CALIBRATION_QUALITIES.get(36), cp.getPhred((char) ('!' + 36), 1));
    } finally {
      assertTrue(FileHelper.deleteAll(dir));
    }
  }

  public void testCurve2() throws IOException {
    final File dir = FileUtils.createTempDir("calMEP", "test");
    try {
      final File cal = FileHelper.resourceToFile("com/rtg/variant/resources/test2.calibration", new File(dir, "test.cal"));
      final Calibrator c = new Calibrator(Calibrator.getCovariateSet(cal), null);
      c.accumulate(cal);
      final CalibratedMachineErrorParams cp = new CalibratedMachineErrorParams(MachineType.ILLUMINA_PE, c, "readgroup1");
      assertEquals(3, cp.getPhred((char) ('!' + 2), 1));
    } finally {
      assertTrue(FileHelper.deleteAll(dir));
    }
  }

  public void testCurveLimit() throws IOException {
    final File dir = FileUtils.createTempDir("calMEP", "test");
    try {
      final File cal = FileHelper.resourceToFile("com/rtg/variant/resources/test.calibration", new File(dir, "test.cal"));
      final Calibrator c = new Calibrator(Calibrator.getCovariateSet(cal), null);
      c.accumulate(cal);
      final CalibratedMachineErrorParams cp = new CalibratedMachineErrorParams(MachineType.ILLUMINA_PE, c, "readgroup1");
      assertEquals((int) TEST_CALIBRATION_QUALITIES.get(2), cp.getPhred((char) ('!' + 2), 1));
      assertEquals((int) TEST_CALIBRATION_QUALITIES.get(36), cp.getPhred((char) ('!' + 36), 1));
      assertEquals((int) TEST_CALIBRATION_QUALITIES.get(63), cp.getPhred((char) ('!' + 63), 1));
      // test that we return the max when asking for higher values
      assertEquals((int) TEST_CALIBRATION_QUALITIES.get(63), cp.getPhred((char) ('!' + 64), 1));
      assertEquals((int) TEST_CALIBRATION_QUALITIES.get(63), cp.getPhred((char) ('!' + 128), 1));

    } finally {
      assertTrue(FileHelper.deleteAll(dir));
    }
  }

  public void testCountsToQ() {
    assertEquals(10, CalibratedMachineErrorParams.countsToEmpiricalQuality(1, 9, 0.01));
    assertEquals(20, CalibratedMachineErrorParams.countsToEmpiricalQuality(1, 99, 0.01));
    assertEquals(30, CalibratedMachineErrorParams.countsToEmpiricalQuality(1, 999, 0.01));
    assertEquals(40, CalibratedMachineErrorParams.countsToEmpiricalQuality(1, 9999, 0.01));
    assertEquals(50, CalibratedMachineErrorParams.countsToEmpiricalQuality(1, 99999, 0.01));
    assertEquals(10, CalibratedMachineErrorParams.countsToEmpiricalQuality(0, 9, 1.0));
    assertEquals(20, CalibratedMachineErrorParams.countsToEmpiricalQuality(0, 99, 1.0));
    assertEquals(30, CalibratedMachineErrorParams.countsToEmpiricalQuality(0, 999, 1.0));
    assertEquals(40, CalibratedMachineErrorParams.countsToEmpiricalQuality(0, 9999, 1.0));
    assertEquals(50, CalibratedMachineErrorParams.countsToEmpiricalQuality(0, 99999, 1.0));
  }

}
