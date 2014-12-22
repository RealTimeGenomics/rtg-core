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
package com.rtg.variant.eval;

import java.io.File;
import java.io.IOException;

import com.rtg.launcher.OutputParams;
import com.rtg.util.test.params.TestParams;
import com.rtg.vcf.VcfUtils;

import junit.framework.TestCase;

/**
 */
public class VcfEvalParamsTest extends TestCase {

  public void testDefaults() {
    final VcfEvalParams params = VcfEvalParams.builder().create();
    assertEquals(VcfUtils.FORMAT_GENOTYPE_QUALITY, params.scoreField());
    assertEquals(RocSortOrder.DESCENDING, params.sortOrder());
    assertNull(params.outputParams());
    assertNull(params.baselineFile());
    assertNull(params.callsFile());
    assertNull(params.templateFile());
    assertNull(params.sampleName());
    assertFalse(params.useAllRecords());
  }

  public void testBuilder() throws IOException {
    VcfEvalParams.VcfEvalParamsBuilder builder = VcfEvalParams.builder();
    builder = builder.name("blah").outputParams(new OutputParams(new File("out"), false, false)).baseLineFile(new File("mutations")).callsFile(new File("calls")).templateFile(new File("template")).maxLength(199);
    builder = builder.scoreField(VcfUtils.QUAL).sortOrder(RocSortOrder.ASCENDING).sampleName("name");
    builder.rtgStats(true);
    final VcfEvalParams params = builder.create();
    assertEquals("blah", params.name());
    assertEquals(VcfUtils.QUAL, params.scoreField());
    assertEquals(RocSortOrder.ASCENDING, params.sortOrder());
    assertEquals("out", params.outputParams().directory().getName());
    assertEquals("mutations", params.baselineFile().getName());
    assertEquals("calls", params.callsFile().getName());
    assertEquals("template", params.templateFile().getName());
    assertEquals("name", params.sampleName());
    assertEquals(199, params.maxLength());
    assertEquals(params.outputParams().directory(), params.directory());
    assertEquals(new File(new File("out"), "bbbb"), params.file("bbbb"));
    assertTrue(params.rtgStats());
  }

  public void testOmnes() {
    new TestParams(VcfEvalParams.class, VcfEvalParams.VcfEvalParamsBuilder.class).check();
  }
}
