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

import java.io.IOException;

import junit.framework.TestCase;

/**
 * Tests the corresponding class.
 */
public class SimpleJoinTest extends TestCase {

  public void test() throws IOException {
    final RegionDataset ds1 = new RegionDataset(new String[0]);
    ds1.regions().add("sequence 0:123-456");
    ds1.regions().add("sequence 1:123-456");
    final IntColumn col = new IntColumn("x");
    col.add(8);
    col.add(8);
    ds1.addColumn(col);
    final RegionDataset ds2 = new RegionDataset(new String[0]);
    ds2.regions().add("sequence 0:123-456");
    ds2.regions().add("sequence 1:123-456");
    final IntColumn col2 = new IntColumn("y");
    col2.add(4);
    col2.add(4);
    ds2.addColumn(col2);
    final SimpleJoin sj = new SimpleJoin(ds1, "prefix");
    sj.process(ds2);
    assertEquals(2, ds2.columns());
    assertEquals("y", ds2.columnName(0));
    assertEquals("prefixx", ds2.columnName(1));
  }
}
