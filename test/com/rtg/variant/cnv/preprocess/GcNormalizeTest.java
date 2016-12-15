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

import com.rtg.reader.ArraySequencesReader;
import com.rtg.util.StringUtils;

import junit.framework.TestCase;

/**
 * Tests the corresponding class.
 */
public class GcNormalizeTest extends TestCase {

  public void test() throws IOException {
    final RegionDataset ds = new RegionDataset(new String[0]);
    ds.regions().add("sequence 0:123-456");
    ds.regions().add("sequence 1:123-456");
    ds.regions().add("sequence 2:123-456");
    final IntColumn col = new IntColumn("x");
    col.add(8);
    col.add(8);
    col.add(8);
    ds.addColumn(col);
    assertEquals(1, ds.columns());
    new AddGc(new ArraySequencesReader(StringUtils.repeat("A", 1000), StringUtils.repeat("G", 1000), StringUtils.repeat("AG", 500))).process(ds);
    assertEquals(3, ds.columns());
    new GcNormalize(0, 1).process(ds);
    assertEquals(4, ds.columns());
    assertEquals("gcnorm(x)", ds.getColumns().get(3).getName());
    assertEquals("0", ds.getColumns().get(1).toString(0));
    assertEquals("334", ds.getColumns().get(1).toString(1));
    assertEquals("167", ds.getColumns().get(1).toString(2));
  }
}
