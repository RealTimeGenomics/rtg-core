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
