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

import com.rtg.variant.MachineErrorParams;
import com.rtg.variant.MachineErrorParamsBuilder;
import com.rtg.variant.realign.AllPaths;
import com.rtg.variant.realign.RealignParams;
import com.rtg.variant.realign.ScoreFastUnderflow;
import com.rtg.variant.realign.ScoreFastUnderflowCG;

import junit.framework.TestCase;

/**
 */
public class ScoreInterfaceMemoTest extends TestCase {

  public void test() throws IOException {
    final ScoreInterfaceMemoInterface memo = new ScoreInterfaceMemo();
    final MachineErrorParamsBuilder builder = MachineErrorParams.builder();
    final RealignParams me = builder.create().realignParams();
    final AllPaths s1 = memo.getScoreInterface(me);
    assertTrue(s1 instanceof ScoreFastUnderflow);
    assertTrue(s1 == memo.getScoreInterface(me));

    final AllPaths s2 = memo.getScoreInterface(builder.create().realignParams());
    assertTrue(s2 instanceof ScoreFastUnderflow);
    assertFalse(s1 == s2);

    final MachineErrorParams complete = builder.errors("complete").create();
    final AllPaths s3 = memo.getScoreInterface(complete.realignParams());
    assertTrue(s3 instanceof ScoreFastUnderflowCG);
    assertTrue(s3 == memo.getScoreInterface(complete.realignParams()));
  }
}
