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
      final DiscordantToolParams def = DiscordantToolParams.builder().mapped(Collections.emptyList()).outputParams(new OutputParams(new File("blah"), false)).genome(new MockReaderParams(1, 1, SequenceMode.BIDIRECTIONAL.codeType())).create();
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
