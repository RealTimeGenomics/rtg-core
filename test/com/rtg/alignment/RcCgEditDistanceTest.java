/*
 * Copyright (c) 2017. Real Time Genomics Limited.
 *
 * Use of this source code is bound by the Real Time Genomics Limited Software Licence Agreement
 * for Academic Non-commercial Research Purposes only.
 *
 * If you did not receive a license accompanying this file, a copy must first be obtained by email
 * from support@realtimegenomics.com.  On downloading, using and/or continuing to use this source
 * code you accept the terms of that license agreement and any amendments to those terms that may
 * be made from time to time by Real Time Genomics Limited.
 */
package com.rtg.alignment;

import java.io.IOException;

import com.rtg.AbstractTest;
import com.rtg.mode.DnaUtils;
import com.rtg.util.InvalidParamsException;
import com.rtg.variant.MachineErrorParamsBuilder;
import com.rtg.variant.realign.RealignParams;
import com.rtg.variant.realign.RealignParamsImplementation;

/**
 */
public class RcCgEditDistanceTest extends AbstractTest {

  public void testAbstract() {
    final CgGotohEditDistanceTest test = new CgGotohEditDistanceTest();
    test.testCGOverlap0Rev();
    test.testCgOverlapRealWorld7();
  }

  private static class TestGotohEditDistance extends CgGotohEditDistance {
    TestGotohEditDistance(int maxShift, RealignParams params) {
      super(maxShift, params, 0);
    }
    @Override
    public int[] calculateEditDistance(byte[] read, int rlen, byte[] template, int zeroBasedStart, int maxScore, int maxShift, boolean cgLeft) {
      return ActionsHelper.build("=====B====X===============NNNNNN==========", 0, 2);
    }
  }

  public void testCGGotohHandling() {
    final RealignParamsImplementation mParams;
    try {
      mParams = new RealignParamsImplementation(new MachineErrorParamsBuilder().errors("cg_test_errors").create());
    } catch (final InvalidParamsException | IOException e) {
      throw new RuntimeException(e);
    }

    final RcEditDistance rced = new RcEditDistance(new TestGotohEditDistance(7, mParams));
    final int[] actions = rced.calculateEditDistance(DnaUtils.encodeString("attcttactanccccaactgagccc     agtaggagta".replaceAll(" ", "")), 35, DnaUtils.encodeString("atactcctacttttgctgggctcagttgggggaagtagaat"), 0, true, 25, 7, true);

    assertEquals(1, ActionsHelper.zeroBasedTemplateStart(actions));

  }
}
