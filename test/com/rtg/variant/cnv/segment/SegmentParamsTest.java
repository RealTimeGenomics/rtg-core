/*
 * Copyright (c) 2016. Real Time Genomics Limited.
 *
 * Use of this source code is bound by the Real Time Genomics Limited Software Licence Agreement
 * for Academic Non-commercial Research Purposes only.
 *
 * If you did not receive a license accompanying this file, a copy must first be obtained by email
 * from support@realtimegenomics.com.  On downloading, using and/or continuing to use this source
 * code you accept the terms of that license agreement and any amendments to those terms that may
 * be made from time to time by Real Time Genomics Limited.
 */

package com.rtg.variant.cnv.segment;

import java.io.File;

import com.rtg.launcher.OutputParams;
import com.rtg.util.test.params.TestParams;

import junit.framework.TestCase;

public class SegmentParamsTest extends TestCase {

  public void testDefaults() {
    final SegmentParams params = SegmentParams.builder().create();
    assertNull(params.outputParams());
  }

  public void testBuilder() {
    SegmentParams.SegmentParamsBuilder builder = SegmentParams.builder();
    builder = builder.name("blah").outputParams(new OutputParams(new File("out"), false)).controlFile(new File("control.bed")).caseFile(new File("case.bed"));
    final SegmentParams params = builder.create();
    assertEquals("blah", params.name());
    assertEquals("out", params.outputParams().directory().getName());
  }

  public void testOmnes() {
    new TestParams(SegmentParams.class, SegmentParams.SegmentParamsBuilder.class).check();
  }
}
