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

package com.rtg.variant.sv.discord.pattern;

import java.io.File;
import java.io.IOException;

import com.rtg.util.TestUtils;
import com.rtg.util.intervals.RegionRestriction;
import com.rtg.util.test.FileHelper;
import com.rtg.util.test.params.TestParams;

import junit.framework.TestCase;

/**
 */
public class BreakpointPatternParamsTest extends TestCase {
  public void testOmnes() {
    new TestParams(BreakpointPatternParams.class, BreakpointPatternParams.Builder.class).check();
  }

  public void testDefaults() {
    final BreakpointPatternParams def = BreakpointPatternParams.builder().create();
    assertEquals(3, def.minDepth());
    assertEquals(500, def.fragmentLength());
    assertEquals(50, def.sameDistance());
    assertEquals(null, def.region());
    assertEquals(null, def.directory());
    assertEquals(null, def.files());
    TestUtils.containsAll(def.toString()
        , "BreakpointPatternParams "
        , " directory=" + null
        , " fragment-length=" + 500
        , " same-distance=" + 50
        , " min-depth=" + 3
    );
  }

  public void testBuilder() {
    final BreakpointPatternParams.Builder builder = BreakpointPatternParams.builder();
    assertEquals(builder, builder.files(null));
    assertEquals(builder, builder.directory(null));
    assertEquals(builder, builder.fragmentLength(1));
    assertEquals(builder, builder.sameDistance(1));
    assertEquals(builder, builder.setMinDepth(1));
    assertEquals(builder, builder.region(new RegionRestriction("chr1:200-399")));
    assertEquals(builder, builder.self());
  }
  public void testFile() throws IOException {
    final File tmp = FileHelper.createTempDirectory();
    try {
      final BreakpointPatternParams params = BreakpointPatternParams.builder().directory(tmp).create();
      assertEquals(new File(tmp, "monkey"), params.file("monkey"));

    } finally {
      FileHelper.deleteAll(tmp);
    }
  }
}
