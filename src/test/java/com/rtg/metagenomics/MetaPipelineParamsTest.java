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
