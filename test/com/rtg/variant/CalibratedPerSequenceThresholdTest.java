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
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.rtg.calibrate.CalibratedPerSequenceExpectedCoverage;
import com.rtg.calibrate.Calibrator;
import com.rtg.util.TestUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.MemoryPrintStream;
import com.rtg.util.test.FileHelper;
import com.rtg.variant.CalibratedPerSequenceThreshold.ThresholdFunction;

import junit.framework.TestCase;

/**
 */
public class CalibratedPerSequenceThresholdTest extends TestCase {


  public void testThresholding() {
    assertEquals(60, ThresholdFunction.SIMPLE_MULT.threshold(30, 2));
    assertEquals(8, ThresholdFunction.SIMPLE_MULT.threshold(4, 2));

    // Under sqrt_mult, dial of 5.5. gives same result as simple_mult
    // with the current default of 2 for 30x coverage, but more robust
    // threshold at low coverage
    assertEquals(60, ThresholdFunction.SQRT_MULT.threshold(30, 5.5));
    assertEquals(15, ThresholdFunction.SQRT_MULT.threshold(4, 5.5));
  }

  public void testCalibratedSimpleMult() throws IOException {
    final MemoryPrintStream mps = new MemoryPrintStream();
    Diagnostic.setLogStream(mps.printStream());
    try {
        final HashMap<String, String> readToSample = new HashMap<>();
        readToSample.put("some_reads", "male");
        final Calibrator c = getCalibrator("com/rtg/variant/resources/calibration_sample.calibration");
        final Map<String, Integer> lengths = new HashMap<>();
        lengths.put("chr1", 10000);
        lengths.put("chr2", 20000);
        lengths.put("chrX", 2000);
        lengths.put("chrY", 1000);
        c.setSequenceLengths(lengths);
        final CalibratedPerSequenceThreshold threshold = new CalibratedPerSequenceThreshold(new CalibratedPerSequenceExpectedCoverage(c, null, readToSample, null), 2, ThresholdFunction.SIMPLE_MULT);
        assertEquals(60, threshold.thresholdSingle("chr1"));
        assertEquals(47, threshold.thresholdSingle("chr2"));
        assertEquals(27, threshold.thresholdSingle("chrX"));
        assertEquals(10, threshold.thresholdSingle("chrY"));

        assertEquals("avgSeqCov*2.0", threshold.toString());
        try {
          threshold.thresholdSingle("I do not exist");
        } catch (final NoTalkbackSlimException e) {
          assertEquals("Unknown sequence: I do not exist", e.getMessage());
        }
    } finally {
      Diagnostic.setLogStream();
    }
  }

  public void testCalibratedSqrtMult() throws IOException {
    final MemoryPrintStream mps = new MemoryPrintStream();
    Diagnostic.setLogStream(mps.printStream());
    try {
        final HashMap<String, String> readToSample = new HashMap<>();
        readToSample.put("some_reads", "male");
        final Calibrator c = getCalibrator("com/rtg/variant/resources/calibration_sample.calibration");
        final Map<String, Integer> lengths = new HashMap<>();
        lengths.put("chr1", 10000);
        lengths.put("chr2", 20000);
        lengths.put("chrX", 2000);
        lengths.put("chrY", 1000);
        final CalibratedPerSequenceThreshold threshold = new CalibratedPerSequenceThreshold(new CalibratedPerSequenceExpectedCoverage(c, lengths, readToSample, null), 2, ThresholdFunction.SQRT_MULT);
        assertEquals(41, threshold.thresholdSingle("chr1"));
        assertEquals(33, threshold.thresholdSingle("chr2"));
        assertEquals(21, threshold.thresholdSingle("chrX"));
        assertEquals(9, threshold.thresholdSingle("chrY"));
        assertEquals("avgSeqCov+(sqrt(avgSeqCov)*2.0)", threshold.toString());
        try {
          threshold.thresholdSingle("I do not exist");
        } catch (final NoTalkbackSlimException e) {
          assertEquals("Unknown sequence: I do not exist", e.getMessage());
        }
    } finally {
      Diagnostic.setLogStream();
    }
  }

  public void testCalibratedMultiSample() throws IOException {
    final MemoryPrintStream mps = new MemoryPrintStream();
    Diagnostic.setLogStream(mps.printStream());
    try {
      final HashMap<String, String> readToSample = new HashMap<>();
      readToSample.put("SAMPLE1", "GENOME1");
      readToSample.put("SAMPLE2", "GENOME2");
      readToSample.put("SAMPLE3", "GENOME3");

      final File top = FileUtils.createTempDir("calibrated", "test");
      try {
        final File g1MatedCal = new File(top, "g1.mated.calibration");
        FileHelper.resourceToFile("com/rtg/variant/resources/genome1_mated.sam.gz.calibration", g1MatedCal);
        final File g1UnmatedCal = new File(top, "g1.unmated.calibration");
        FileHelper.resourceToFile("com/rtg/variant/resources/genome1_unmated.sam.gz.calibration", g1UnmatedCal);
        final File g2MatedCal = new File(top, "g2.mated.calibration");
        FileHelper.resourceToFile("com/rtg/variant/resources/genome2_mated.sam.gz.calibration", g2MatedCal);
        final File g2UnmatedCal = new File(top, "g2.unmated.calibration");
        FileHelper.resourceToFile("com/rtg/variant/resources/genome2_unmated.sam.gz.calibration", g2UnmatedCal);
        final File g3MatedCal = new File(top, "g3.mated.calibration");
        FileHelper.resourceToFile("com/rtg/variant/resources/genome3_mated.sam.gz.calibration", g3MatedCal);
        final File g3UnmatedCal = new File(top, "g3.unmated.calibration");
        FileHelper.resourceToFile("com/rtg/variant/resources/genome3_unmated.sam.gz.calibration", g3UnmatedCal);
        final List<File> f = new ArrayList<>();
        f.add(g1MatedCal);
        f.add(g1UnmatedCal);
        f.add(g2MatedCal);
        f.add(g2UnmatedCal);
        f.add(g3MatedCal);
        f.add(g3UnmatedCal);
        final Calibrator c = Calibrator.initCalibrator(f);
        final Map<String, Integer> lengths = new HashMap<>();
        lengths.put("simulatedSequence1", 10000);
        lengths.put("simulatedSequence2", 10000);
        lengths.put("simulatedSequence3", 10000);
        c.setSequenceLengths(lengths);
        final CalibratedPerSequenceThreshold threshold = new CalibratedPerSequenceThreshold(new CalibratedPerSequenceExpectedCoverage(c, lengths, readToSample, null), 2, ThresholdFunction.SIMPLE_MULT);
        assertEquals(29, threshold.thresholdSingle("simulatedSequence1"));
        assertEquals(31, threshold.thresholdSingle("simulatedSequence2"));
        assertEquals(30, threshold.thresholdSingle("simulatedSequence3"));
        assertEquals(58, threshold.thresholdTotal("simulatedSequence1"));
        assertEquals(61, threshold.thresholdTotal("simulatedSequence2"));
        assertEquals(60, threshold.thresholdTotal("simulatedSequence3"));
        assertEquals("avgSeqCov*2.0", threshold.toString());

      } finally {
        assertTrue(FileHelper.deleteAll(top));
      }
    } finally {
      Diagnostic.setLogStream();
    }
  }

  public static Calibrator getCalibrator(String resource) throws IOException {

    final File calibration = File.createTempFile("CalibratedPerSequenceThresholdTest", ".calibration");
    try {
      FileHelper.resourceToFile(resource, calibration);
      final Calibrator c = new Calibrator(Calibrator.getCovariateSet(calibration), null);
      c.accumulate(calibration);
      return c;
    } finally {
      FileHelper.deleteAll(calibration);
    }
  }

  public void testMissingCovariate() throws IOException {
    final MemoryPrintStream mps = new MemoryPrintStream();
    Diagnostic.setLogStream(mps.printStream());
    try {
      final HashMap<String, String> readToSample = new HashMap<>();
      readToSample.put("some_reads", "male");

      final Calibrator c = getCalibrator("com/rtg/variant/resources/calibration_global_sample.calibration");
      final Map<String, Integer> lengths = new HashMap<>();
      lengths.put("chr1", 10000);
      lengths.put("chr2", 20000);
      lengths.put("chrX", 2000);
      lengths.put("chrY", 1000);
      c.setSequenceLengths(lengths);
      final CalibratedPerSequenceThreshold threshold = new CalibratedPerSequenceThreshold(new CalibratedPerSequenceExpectedCoverage(c, lengths, readToSample, null), 2, ThresholdFunction.SIMPLE_MULT);
      assertEquals(20, threshold.thresholdSingle("chr1"));
      assertEquals(20, threshold.thresholdSingle("chr2"));
      assertEquals(20, threshold.thresholdSingle("chrX"));
      assertEquals(20, threshold.thresholdSingle("chrY"));
      assertEquals("avgSeqCov*2.0", threshold.toString());
      try {
        threshold.thresholdSingle("I do not exist");
      } catch (final NoTalkbackSlimException e) {
        assertEquals("Unknown sequence: I do not exist", e.getMessage());
      }
    } finally {
      Diagnostic.setLogStream();
    }
    TestUtils.containsAll(mps.toString(),
      "Average coverage for sample male is 10.09"
    );
  }

}
