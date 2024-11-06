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
package com.rtg.variant.cnv;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

import com.rtg.launcher.OutputParams;
import com.rtg.util.TestUtils;
import com.rtg.util.io.FileUtils;
import com.rtg.util.test.params.TestParams;
import com.rtg.variant.cnv.CnvProductParams.CnvProductParamsBuilder;

import junit.framework.TestCase;

/**
 * Test class
 */
public class CnvProductParamsTest extends TestCase {

  /** Because we want it to be hard to be agile! yay */
  public void testDefaults() {
    final CnvProductParams.CnvProductParamsBuilder builder = CnvProductParams.builder();
    assertNotNull(builder);
    final CnvProductParams params = builder.create();
    assertNotNull(params);
    assertEquals(100, params.bucketSize());
    assertNull(params.directory());

    try {
      params.file("file");
      fail();
    } catch (final RuntimeException e) {
      //ignore
    }
    assertEquals(null, params.outputParams());
    assertEquals(true, params.filterStartPositions());
    assertEquals(null, params.mappedBase());
    assertEquals(null, params.mappedTarget());
    assertEquals(-1, params.filterParams().maxAlignmentCount());
    assertEquals(null, params.filterParams().maxMatedAlignmentScore());
    assertEquals(null, params.filterParams().maxUnmatedAlignmentScore());
    assertEquals(3.0, params.divisionFactor());
    assertEquals(3.0, params.multiplicationFactor());
    TestUtils.containsAll(params.toString(), "CnvProductParams",
      "bucketSize=100",
      "baseLineInput=null",
      "targetInput=null",
      "filterStartPositions=" + true,
      "divisionFactor=" + 3.0,
      "multiplicationFactor=" + 3.0,
      "null",
      "maxMatedAlignmentScore=null",
      "maxUnmatedAlignmentScore=null",
      "maxAlignmentCount=-1"
      );
  }

  public void testActual() throws IOException {
    final File tempFile = FileUtils.createTempDir("cnvproductparams", "test");
    assertTrue(tempFile.delete());
    final CnvProductParams params = CnvProductParams.builder()
      .mappedBase(new ArrayList<>())
      .mappedTarget(new ArrayList<>())
      .outputParams(new OutputParams(tempFile, false))
      .divisionFactor(5.0)
      .multiplicationFactor(7.0)
      .bucketSize(2)
      .filterStartPositions(false)
      .magicConstant(1.3)
      .extraPenaltyOff(true)
      .threads(9)
      .create()
      ;
    TestUtils.containsAll(params.toString(), "CnvProductParams",
      "bucketSize=2",
      "baseLineInput=[]",
      "targetInput=[]",
      "filterStartPositions=" + false,
      "divisionFactor=" + 5.0,
      "multiplicationFactor=" + 7.0,
      "OutputParams",
      "directory=" + tempFile.getAbsolutePath(),
      "maxMatedAlignmentScore=null",
      "maxUnmatedAlignmentScore=null",
      "maxAlignmentCount=-1",
      "zip=" + false);


    assertEquals(1.3, params.magicConstant(), 1e-8);
    assertTrue(params.extraPenaltyOff());
    assertEquals(9, params.threads());
    assertEquals(tempFile, params.directory());
    assertEquals(new File(tempFile, "file").getPath(), params.file("file").getPath());
  }

  public void testOmnes() {
    new TestParams(CnvProductParams.class, CnvProductParamsBuilder.class).check();
  }

}
