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
package com.rtg.variant.avr;

import java.io.File;

import junit.framework.TestCase;

/**
 */
public class VcfDatasetTest extends TestCase {

  public void testAnnotate() {
    final File f = new File("boo");
    final VcfDataset d = new VcfDataset(f, 5, true, false, 1.5);
    assertEquals(f, d.getVcfFile());
    assertEquals(5, d.getSampleNum());
    assertTrue(d.isPositive());
    assertEquals(1.5, d.getInstanceWeight());
  }
}
