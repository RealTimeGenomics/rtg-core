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
package com.rtg.variant.realign;

import com.rtg.util.integrity.Exam;
import com.rtg.util.machine.MachineType;
import com.rtg.variant.AbstractMachineErrorParams;
import com.rtg.variant.MachineErrorParams;

import junit.framework.TestCase;

/**
 */
public class RealignParamsImplementationTest extends TestCase {

  public void testWeightedDecay() {
    assertEquals(-1.62985, RealignParamsImplementation.weightedDecay(new double[] {0.0, 0.8, 0.15, 0.03, 0.01, 0.01}), 0.00001);
    assertEquals(-1.38629, RealignParamsImplementation.weightedDecay(new double[] {0.0, 0.8, 0.2}), 0.00001);
    assertEquals(-1.60944, RealignParamsImplementation.weightedDecay(new double[] {0.0, 1.0}), 0.00001);
    assertEquals(-1.60944, RealignParamsImplementation.weightedDecay(new double[] {0.0, 0.0, 1.0}), 0.00001);
  }

  public void testMean() {
    assertEquals(1.0, RealignParamsImplementation.mean(new double[] {0.0, 1.0}), 0.00001);
    assertEquals(1.5, RealignParamsImplementation.mean(new double[] {0.0, 0.5, 0.5}), 0.00001);
  }

  public void test() {
    final AbstractMachineErrorParams params = MachineErrorParams.builder()
    .errorInsEventRate(0.1)
    .errorInsDistribution(new double[] {0.0, 0.8, 0.15, 0.03, 0.01, 0.01})
    .errorDelEventRate(0.2)
    .errorDelDistribution(new double[] {0.0, 0.8, 0.2})
    .errorMnpEventRate(0.3)
    .errorMnpDistribution(new double[] {0.5, 0.5})
    .create();
    final RealignParams rp = new RealignParamsImplementation(params);
    Exam.integrity(rp);
    assertEquals(Math.log(0.1), rp.insertOpenLn(), 0.00001);
    assertEquals(-1.62985, rp.insertExtendLn(), 0.00001);
    assertEquals(Math.log(0.2), rp.deleteOpenLn(), 0.00001);
    assertEquals(-1.38629, rp.deleteExtendLn(), 0.00001);
    assertEquals(Math.log(0.3 * 1.5), rp.misMatchLn(), 0.00001);
    assertEquals(Math.log(1.0 - 0.3 * 1.5), rp.matchLn(), 0.00001);
  }

  public void testCG() {
    final AbstractMachineErrorParams params = MachineErrorParams.builder().create();
    final RealignParams rp = new RealignParamsImplementation(params);
    assertEquals(false, rp.machineType() == MachineType.COMPLETE_GENOMICS);
    assertEquals(-3, rp.gapStart(RealignParamsImplementation.CG_OVERLAP));
    assertEquals(-1, rp.gapEnd(RealignParamsImplementation.CG_OVERLAP));
    assertEquals(0.08, Math.exp(rp.gapFreqLn(RealignParamsImplementation.CG_OVERLAP, -3)), 1E-16);
    assertEquals(0.84, Math.exp(rp.gapFreqLn(RealignParamsImplementation.CG_OVERLAP, -2)), 1E-16);
    assertEquals(0.08, Math.exp(rp.gapFreqLn(RealignParamsImplementation.CG_OVERLAP, -1)), 1E-16);
    //  // NOTE that the small and large gap were back to front here, until about revision 30470
    //    {Math.log(0.08), Math.log(0.84), Math.log(0.08)},    // for gap of -3, -2, -1
    //    {Math.log(0.27), Math.log(0.635), Math.log(0.095)},  // for gap of  0,  1,  2
    //    {Math.log(0.90), Math.log(0.07), Math.log(0.03)},    // for gap of  5,  6,  7
    assertEquals(0, rp.gapStart(RealignParamsImplementation.CG_SMALL_GAP));
    assertEquals(2, rp.gapEnd(RealignParamsImplementation.CG_SMALL_GAP));
    assertEquals(0.90, Math.exp(rp.gapFreqLn(RealignParamsImplementation.CG_SMALL_GAP, 0)), 1E-16);
    assertEquals(0.07, Math.exp(rp.gapFreqLn(RealignParamsImplementation.CG_SMALL_GAP, 1)), 1E-16);
    assertEquals(0.03, Math.exp(rp.gapFreqLn(RealignParamsImplementation.CG_SMALL_GAP, 2)), 1E-16);

    assertEquals(5, rp.gapStart(RealignParamsImplementation.CG_LARGE_GAP));
    assertEquals(7, rp.gapEnd(RealignParamsImplementation.CG_LARGE_GAP));
    assertEquals(0.27, Math.exp(rp.gapFreqLn(RealignParamsImplementation.CG_LARGE_GAP, 5)), 1E-16);
    assertEquals(0.64, Math.exp(rp.gapFreqLn(RealignParamsImplementation.CG_LARGE_GAP, 6)), 1E-16);
    assertEquals(0.09, Math.exp(rp.gapFreqLn(RealignParamsImplementation.CG_LARGE_GAP, 7)), 1E-16);
  }
}
