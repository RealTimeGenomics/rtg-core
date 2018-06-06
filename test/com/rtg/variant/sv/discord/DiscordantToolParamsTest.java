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

package com.rtg.variant.sv.discord;

import java.io.File;
import java.util.Collections;

import com.rtg.launcher.MockReaderParams;
import com.rtg.launcher.OutputParams;
import com.rtg.mode.SequenceMode;
import com.rtg.util.TestUtils;
import com.rtg.util.test.params.TestParams;
import com.rtg.variant.sv.discord.DiscordantToolParams.DiscordantToolParamsBuilder;

import junit.framework.TestCase;

/**
 */
public class DiscordantToolParamsTest extends TestCase {

  public void testOmnes() {
    new TestParams(DiscordantToolParams.class, DiscordantToolParamsBuilder.class).check();
  }

  public void testDefaults() {
      final DiscordantToolParams def = DiscordantToolParams.builder().mapped(Collections.emptyList()).outputParams(new OutputParams(new File("blah"), false, false)).genome(new MockReaderParams(1, 1, SequenceMode.BIDIRECTIONAL)).create();
      assertFalse(def.bedOutput());
      assertFalse(def.debugOutput());
      assertFalse(def.intersectionOnly());
      assertTrue(def.outputTabixIndex());
      assertEquals(0, def.minBreakpointDepth());
      TestUtils.containsAll(def.toString()
          , "DiscordantToolParams "
          , " bed-output=" + false
          , " debug-output=" + false
          , " intersection-only=" + false
          , " output-tabix-index=" + true
          , " min-breakpoint-depth=" + 0
          );
  }

  public void testBuilder() {
    final DiscordantToolParamsBuilder builder = DiscordantToolParams.builder();
    assertEquals(builder, builder.bedOutput(true));
    assertEquals(builder, builder.debugOutput(true));
    assertEquals(builder, builder.intersectionOnly(true));
    assertEquals(builder, builder.outputTabixIndex(false));
    assertEquals(builder, builder.minBreakpointDepth(1));
  }
}
