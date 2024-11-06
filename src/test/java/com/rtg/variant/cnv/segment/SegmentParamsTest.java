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
    builder = builder.name("blah")
      .outputParams(new OutputParams(new File("out"), false))
      .controlFile(new File("control.bed"))
      .caseFile(new File("case.bed"))
      .panelFile(new File("panel.bed"))
      .summaryRegionsFile(new File("genes.bed"))
      .sampleName("foo")
      .coverageColumnName("bar")
      .panelCoverageColumnName("baz")
      .precomputedColumn(53)
      .minBins(54)
      .gcBins(55)
      .minSegments(56)
      .maxSegments(57)
      .aleph(0.1234)
      .alpha(0.1235)
      .beta(0.1236)
      .minLogR(0.1237)
      .minNormControlCoverage(0.1238)
      .minControlCoverage(0.1239)
      .minCaseCoverage(0.1240)
      .absorbSingletons(true)
      .graphviz(true)
    ;
    final SegmentParams params = builder.create();
    assertEquals("blah", params.name());
    assertEquals("out", params.outputParams().directory().getName());
    assertEquals("control.bed", params.controlFile().getName());
    assertEquals("case.bed", params.caseFile().getName());
    assertEquals("panel.bed", params.panelFile().getName());
    assertEquals("genes.bed", params.summaryRegionsFile().getName());
    assertEquals("foo", params.sampleName());
    assertEquals("bar", params.coverageColumnName());
    assertEquals("baz", params.panelCoverageColumnName());
    assertEquals(53, params.precomputedColumn());
    assertEquals(54, params.minBins());
    assertEquals(55, params.gcBins());
    assertEquals(56, params.minSegments());
    assertEquals(57, params.maxSegments());
    assertEquals(0.1234, params.aleph());
    assertEquals(0.1235, params.alpha());
    assertEquals(0.1236, params.beta());
    assertEquals(0.1237, params.minLogR());
    assertEquals(0.1238, params.minNormControlCoverage());
    assertEquals(0.1239, params.minControlCoverage());
    assertEquals(0.1240, params.minCaseCoverage());
    assertTrue(params.absorbSingletons());
    assertTrue(params.graphviz());
  }

  public void testOmnes() {
    new TestParams(SegmentParams.class, SegmentParams.SegmentParamsBuilder.class).check();
  }
}
