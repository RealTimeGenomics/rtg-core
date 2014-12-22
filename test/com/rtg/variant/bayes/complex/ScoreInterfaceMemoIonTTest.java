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

import java.io.IOException;
import java.io.StringReader;

import com.rtg.util.InvalidParamsException;
import com.rtg.variant.MachineErrorParams;
import com.rtg.variant.realign.HomoPolymerParams;
import com.rtg.variant.realign.RealignParams;
import com.rtg.variant.realign.RealignParamsImplementation;
import com.rtg.variant.util.arithmetic.LogApproximatePossibility;
import com.rtg.variant.util.arithmetic.SimplePossibility;

import junit.framework.TestCase;

/**
 * Test class
 */
public class ScoreInterfaceMemoIonTTest extends TestCase {

  public void test() throws IOException, InvalidParamsException {
    final ScoreInterfaceMemoIonT memo = new ScoreInterfaceMemoIonT(new HomoPolymerParams(SimplePossibility.SINGLETON, 2, 2, new StringReader("")), new HomoPolymerParams(LogApproximatePossibility.SINGLETON, 2, 2, new StringReader("")));
    final RealignParams p = new RealignParamsImplementation(MachineErrorParams.builder().create());
    assertTrue(memo.getScoreInterface(p, false) == memo.getScoreInterface(p, false));
    try {
      memo.getScoreInterface(p, true);
      fail();
    } catch (final UnsupportedOperationException e) {

    }
  }
}
