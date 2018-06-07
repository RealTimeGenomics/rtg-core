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
      .outputParams(new OutputParams(tempFile, false, false))
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
      "progress=" + false,
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
