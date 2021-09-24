/*
 * Copyright (c) 2018. Real Time Genomics Limited.
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

/**
 * Tests the corresponding class.
 */
public class UnionJoinTest extends SimpleJoinTest {

  @Override
  protected DatasetProcessor getJoiner(RegionDataset d) {
    return new UnionJoin(d, "prefix");
  }

  public void test2() throws IOException {
    final RegionDataset ds1 = new RegionDataset(new String[0]);
    ds1.regions().add("sequence 0:123-456");
    ds1.regions().add("sequence 1:123-456");
    ds1.regions().add("sequence 2:123-456");
    final IntColumn col = new IntColumn("x");
    col.add(8);
    col.add(10);
    col.add(12);
    ds1.addColumn(col);
    final RegionDataset ds2 = new RegionDataset(new String[0]);
    ds2.regions().add("sequence 0:123-456");
    ds2.regions().add("sequence 1:123-456");
    ds2.regions().add("sequence 1:789-101112");
    final IntColumn col2 = new IntColumn("y");
    col2.add(2);
    col2.add(4);
    col2.add(6);
    ds2.addColumn(col2);
    final DatasetProcessor sj = getJoiner(ds1);
    sj.process(ds2);
    assertEquals(2, ds2.columns());
    assertEquals("y", ds2.columnName(0));
    assertEquals("prefixx", ds2.columnName(1));
    assertEquals(4, ds2.size());
    assertEquals(6.0, ((NumericColumn) ds2.column(0)).get(2));
    assertTrue(Double.isNaN(((NumericColumn) ds2.column(0)).get(3)));
    assertTrue(Double.isNaN(((NumericColumn) ds2.column(1)).get(2)));
    assertEquals(12.0, ((NumericColumn) ds2.column(1)).get(3));
  }

  public void test3() throws IOException {
    final RegionDataset ds1 = new RegionDataset(new String[0]);
    final RegionDataset ds2 = new RegionDataset(new String[0]);
    ds1.regions().add("sequence 0:123-456");
    ds1.regions().add("sequence 1:123-456");
    ds1.regions().add("sequence 1:789-101112");
    final DatasetProcessor sj = getJoiner(ds1);
    sj.process(ds2);
    assertEquals(3, ds2.size());
  }

  public void test4() throws IOException {
    final RegionDataset ds1 = new RegionDataset(new String[0]);
    ds1.regions().add("sequence 0:123-456");
    ds1.regions().add("sequence 2:123-456");
    ds1.regions().add("sequence 3:789-101112");
    ds1.regions().add("sequence 6:789-101112");
    final DatasetProcessor sj = getJoiner(ds1);

    final RegionDataset ds2 = new RegionDataset(new String[0]);
    ds2.regions().add("sequence 0:123-456");
    ds2.regions().add("sequence 5:123-456");
    ds2.regions().add("sequence 6:789-101112");

    sj.process(ds2);
    assertEquals(5, ds2.size());
    assertEquals("sequence 0", ds2.regions().get(0).getSequenceName());
    assertEquals("sequence 5", ds2.regions().get(1).getSequenceName());
    assertEquals("sequence 2", ds2.regions().get(2).getSequenceName());
    assertEquals("sequence 3", ds2.regions().get(3).getSequenceName());
    assertEquals("sequence 6", ds2.regions().get(4).getSequenceName());
  }
}
