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

package com.rtg.metagenomics;

import java.io.File;

import com.rtg.metagenomics.MetaPipelineParams.MetaPipelineParamsBuilder;
import com.rtg.metagenomics.MetagenomicsWrapperCli.Platform;

import junit.framework.TestCase;

/**
 */
public class MetaPipelineParamsTest extends TestCase {

  public void testDefaults() {
    final MetaPipelineParams params = MetaPipelineParams.builder().create();
    assertNull(params.inputFile());
    assertNull(params.inputLeft());
    assertNull(params.inputRight());
    assertNull(params.filterSdf());
    assertNull(params.speciesSdf());
    assertNull(params.proteinSdf());
    assertEquals(Platform.ILLUMINA, params.inputPlatform());
  }

  public void testBuilder() {
    final MetaPipelineParamsBuilder builder = MetaPipelineParams.builder();
    assertEquals(builder, builder.inputFile(new File("input")));
    assertEquals(builder, builder.inputLeft(new File("left")));
    assertEquals(builder, builder.inputRight(new File("right")));
    assertEquals(builder, builder.filterSdf(new File("filter")));
    assertEquals(builder, builder.speciesSdf(new File("species")));
    assertEquals(builder, builder.proteinSdf(new File("protein")));
    assertEquals(builder, builder.inputPlatform(Platform.IONTORRENT));
    final MetaPipelineParams params = builder.create();
    assertEquals("input", params.inputFile().getName());
    assertEquals("left", params.inputLeft().getName());
    assertEquals("right", params.inputRight().getName());
    assertEquals("filter", params.filterSdf().getName());
    assertEquals("species", params.speciesSdf().getName());
    assertEquals("protein", params.proteinSdf().getName());
    assertEquals(Platform.IONTORRENT, params.inputPlatform());
  }
}
