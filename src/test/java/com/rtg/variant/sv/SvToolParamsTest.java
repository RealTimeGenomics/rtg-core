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
package com.rtg.variant.sv;




import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.util.Collections;

import com.rtg.launcher.MockReaderParams;
import com.rtg.launcher.OutputParams;
import com.rtg.mode.SequenceMode;
import com.rtg.util.TestUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.FileUtils;
import com.rtg.util.test.FileHelper;
import com.rtg.util.test.params.TestParams;
import com.rtg.variant.sv.SvToolParams.SvToolParamsBuilder;

import junit.framework.TestCase;

/**
 */
public class SvToolParamsTest extends TestCase {

  public void testOmnes() {
    new TestParams(SvToolParams.class, SvToolParamsBuilder.class).check();
  }

  public void testDefaults() throws IOException {
    Diagnostic.setLogStream();
    final File tempDir = FileUtils.createTempDir("svtoolparams", "test");
    try {
      final SvToolParams def = SvToolParams.builder().mapped(Collections.emptyList()).outputParams(new OutputParams(tempDir, false)).genome(new MockReaderParams(1, 1, SequenceMode.BIDIRECTIONAL.codeType())).create();
      assertEquals(def.binSize(), 0);
      assertEquals(def.stepSize(), 0);
      assertEquals(def.fineStepSize(), 0);
      assertFalse(def.heterozygous());
      assertFalse(def.outputSimple());
      assertNull(def.correctionsFile());
      try (OutputStream bayStream = def.bayesianStream()) {
        bayStream.write("blah blah".getBytes());
      }
      final File bayFile = new File(tempDir, "sv_bayesian.tsv");
      assertTrue(bayFile.isFile());
      assertEquals("blah blah", FileUtils.fileToString(bayFile));
      TestUtils.containsAll(def.toString()
          , "SvToolParams "
          , " output-simple=" + false
          , " heterozygous=" + false
          , " bin-size=" + 0
          , " step-size=" + 0
          , " fine-step-size=" + 0
          , " correctionsFile=" + null
          );
    } finally {
      assertTrue(FileHelper.deleteAll(tempDir));
    }
  }

  public void testBuilder() {
    final SvToolParamsBuilder builder = SvToolParams.builder();
    assertEquals(builder, builder.binSize(1));
    assertEquals(builder, builder.stepSize(1));
    assertEquals(builder, builder.fineStepSize(1));
    assertEquals(builder, builder.heterozygous(true));
    assertEquals(builder, builder.outputSimple(true));
    assertEquals(builder, builder.correctionsFile(null));
  }
}
