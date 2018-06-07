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
package com.rtg.position.output;

import java.io.File;

import com.rtg.util.TestUtils;
import com.rtg.util.Utils;

import junit.framework.TestCase;

/**
 */
public class PositionOutputParamsTest extends TestCase {

  PositionOutputParams getParams(final File outdir, final OutputFormatType format, final Integer maxGap, final boolean zip, final int topN) {
    final PositionDistributionParams distr = new PositionDistributionParams(0.001, 0.009, maxGap, 0);
    return new PositionOutputParams(outdir, format, maxGap == null ? null : distr, null, zip, topN);
  }

  public void testEquals() {
    final File fa = new File("a");
    final File fb = new File("b");

    final PositionOutputParams a1 = getParams(fa, OutputFormatType.SEGMENT, null, false, 1);
    a1.integrity();
    assertEquals(null, a1.maxGap());
    final PositionOutputParams a2 = getParams(fa, OutputFormatType.SEGMENT, null, false, 1);
    a2.integrity();
    final PositionOutputParams e = getParams(fb, OutputFormatType.SEGMENT, 42, false, 1);
    e.integrity();
    assertEquals(Integer.valueOf(42), e.maxGap());
    final PositionOutputParams f = getParams(fb, OutputFormatType.SEGMENT, 43, false, 1);
    f.integrity();
    final PositionOutputParams j = getParams(fb, OutputFormatType.SEGMENT, 43, true, 1);
    j.integrity();
    TestUtils.equalsHashTest(new PositionOutputParams[][] {{a1, a2}, {e}, {f}, {j}});
  }

  public void testDistr() {
    final File fa = new File("a");
    final PositionDistributionParams distr = new PositionDistributionParams(0.001, 0.009, 42, 0);
    final PositionOutputParams pp = new PositionOutputParams(fa, OutputFormatType.SEGMENT, distr, 0.0, false, 0);
    assertEquals(distr, pp.distribution());
  }

  public void test() {
    final File fa = new File("a");
    final PositionOutputParams a1 = getParams(fa, OutputFormatType.SEGMENT, null, false, 1);
    assertEquals("outputDir=a format=SEGMENT zip=" + Boolean.FALSE.toString() + " score threshold=" + Utils.realFormat(null) + " topN=1"  + " distribution={null}", a1.toString());
    assertEquals(fa, a1.directory());
    assertEquals(OutputFormatType.SEGMENT, a1.format());
    assertEquals(1, a1.topN());
  }
}

