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

import com.rtg.reader.ArraySequencesReader;
import com.rtg.util.StringUtils;

import junit.framework.TestCase;

/**
 * Tests the corresponding class.
 */
public class AddGcTest extends TestCase {

  public void test() throws IOException {
    final RegionDataset ds = new RegionDataset(new String[0]);
    ds.regions().add("sequence 0:123-456");
    ds.regions().add("sequence 1:123-456");
    ds.regions().add("sequence 2:123-456");
    new AddGc(new ArraySequencesReader(StringUtils.repeat("A", 1000), StringUtils.repeat("G", 1000), StringUtils.repeat("AG", 500))).process(ds);
    assertEquals(2, ds.columns());
    assertEquals("gc_content_rel", ds.getColumns().get(1).getName());
    assertEquals("0.00000", ds.getColumns().get(1).toString(0));
    assertEquals("1.00000", ds.getColumns().get(1).toString(1));
    assertEquals("0.500000", ds.getColumns().get(1).toString(2));
  }
}
