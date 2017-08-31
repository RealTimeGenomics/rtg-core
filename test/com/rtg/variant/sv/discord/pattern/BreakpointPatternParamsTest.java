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
