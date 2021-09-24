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

package com.rtg.variant.cnv.preprocess;

import junit.framework.TestCase;

/**
 * Tests the corresponding class.
 */
public class AddLengthTest extends TestCase {

  public void test() {
    final RegionDataset ds = new RegionDataset(new String[0]);
    ds.regions().add("1:42+1");
    ds.regions().add("2:42+10");
    new AddLength().process(ds);
    assertEquals(1, ds.columns());
    assertEquals("length", ds.getColumns().get(0).getName());
    assertEquals("1", ds.getColumns().get(0).toString(0));
    assertEquals("10", ds.getColumns().get(0).toString(1));
  }
}
