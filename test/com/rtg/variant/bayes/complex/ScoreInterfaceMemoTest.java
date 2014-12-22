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
package com.rtg.variant.bayes.complex;

import com.rtg.variant.AbstractMachineErrorParams;
import com.rtg.variant.MachineErrorParams;
import com.rtg.variant.realign.AllPaths;
import com.rtg.variant.realign.RealignParams;
import com.rtg.variant.realign.RealignParamsImplementation;
import com.rtg.variant.realign.ScoreFastUnderflow;
import com.rtg.variant.realign.ScoreFastUnderflowCG;

import junit.framework.TestCase;

/**
 */
public class ScoreInterfaceMemoTest extends TestCase {

  public void test() {
    final ScoreInterfaceMemoInterface memo = new ScoreInterfaceMemo();
    final AbstractMachineErrorParams params = MachineErrorParams.builder()
      .errorInsEventRate(0.1)
      .errorInsDistribution(new double[] {0.0, 0.8, 0.15, 0.03, 0.01, 0.01})
      .errorDelEventRate(0.2)
      .errorDelDistribution(new double[] {0.0, 0.8, 0.2})
      .errorMnpEventRate(0.3)
      .errorMnpDistribution(new double[] {0.5, 0.5})
      .create();
    final RealignParams p1 = new RealignParamsImplementation(params);
    final AllPaths s1 = memo.getScoreInterface(p1, false);
    assertTrue(s1 instanceof ScoreFastUnderflow);
    assertTrue(s1 == memo.getScoreInterface(p1, false));
    final RealignParams p2 = new RealignParamsImplementation(params);
    final AllPaths s2 = memo.getScoreInterface(p2, false);
    assertTrue(s2 instanceof ScoreFastUnderflow);
    assertFalse(s1 == s2);
    final AllPaths s3 = memo.getScoreInterface(p1, true);
    assertTrue(s3 instanceof ScoreFastUnderflowCG);
    assertTrue(s3 == memo.getScoreInterface(p1, true));
  }
}
