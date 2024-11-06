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

