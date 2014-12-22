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

package com.rtg.util.machine;

import com.rtg.util.TestUtils;

import junit.framework.TestCase;

/**
 */
public class MachineTypeTest extends TestCase {

  public void test() {
    TestUtils.testPseudoEnum(MachineType.class, "[illumina_se, illumina_pe, complete_genomics, 454_pe, 454_se, iontorrent]");
  }

  public void testCompatible() {
    assertTrue(MachineType.ILLUMINA_SE.compatiblePlatform("illumina"));
    assertTrue(MachineType.ILLUMINA_SE.compatiblePlatform("Illumina"));

    assertTrue(MachineType.ILLUMINA_PE.compatiblePlatform("illumina"));
    assertTrue(MachineType.ILLUMINA_PE.compatiblePlatform("Illumina"));

    assertTrue(MachineType.COMPLETE_GENOMICS.compatiblePlatform("Complete"));

    assertTrue(MachineType.FOURFIVEFOUR_PE.compatiblePlatform("lS454"));
    assertTrue(MachineType.FOURFIVEFOUR_SE.compatiblePlatform("Ls454"));

    assertTrue(MachineType.IONTORRENT.compatiblePlatform("iontorrent"));
  }
}
