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
package com.rtg.variant.cnv.preprocess;

import java.io.IOException;

/**
 * Tests the corresponding class.
 */
public class IntersectJoinTest extends SimpleJoinTest {

  @Override
  protected DatasetProcessor getJoiner(RegionDataset d) {
    return new IntersectJoin(d, "prefix");
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
    assertEquals(2, ds2.size());
  }
}
