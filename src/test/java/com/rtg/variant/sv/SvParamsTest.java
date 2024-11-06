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

import com.rtg.launcher.OutputParams;
import com.rtg.sam.SamFilterParams;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.FileUtils;
import com.rtg.util.test.FileHelper;
import com.rtg.util.test.params.TestParams;
import com.rtg.variant.sv.SvParams.SvParamsBuilder;
import com.rtg.variant.sv.SvParamsTest.MockSvParams.MockSvParamsBuilder;

import junit.framework.TestCase;

/**
 */
public class SvParamsTest extends TestCase {

  static final class MockSvParams extends SvParams {

    static MockSvParamsBuilder builder() {
      return new MockSvParamsBuilder();
    }

    static final class MockSvParamsBuilder extends SvParamsBuilder<MockSvParamsBuilder> {
      @Override
      protected MockSvParamsBuilder self() {
        return this;
      }

      public MockSvParams create() {
        return new MockSvParams(this);
      }
    }

    MockSvParams(MockSvParamsBuilder builder) {
      super(builder);
    }
  }

  public void testOmnes() {
    new TestParams(SvParams.class, SvParamsBuilder.class).check();
  }

  public void testDefaultSvParams() throws IOException {
    Diagnostic.setLogStream();
    final File tempDir = FileUtils.createTempDir("svparams", "test");
    try {
      final MockSvParams def = MockSvParams.builder().outputParams(new OutputParams(tempDir, true)).create();

      assertNotNull(def.outputParams());
      assertNotNull(def.filterParams());
      assertNull(def.readGroupLabels());
      assertNull(def.readGroupStatistics());
      assertTrue(def.blockCompressed());
      assertEquals(tempDir, def.directory());
      assertEquals(new File(tempDir, "blah.txt.gz"), def.outFile("blah.txt"));
      assertEquals(new File(tempDir, "blah.txt"), def.file("blah.txt"));
      try (OutputStream str = def.outStream("df.txt")) {
        str.write("blah blah".getBytes());
      }
      final File file = new File(tempDir, "df.txt.gz");
      assertTrue(file.exists());
      assertEquals("blah blah", FileHelper.gzFileToString(file));
    } finally {
      assertTrue(FileHelper.deleteAll(tempDir));
    }
  }

  public void testSvParamsBuilder() {
    final MockSvParamsBuilder builder = MockSvParams.builder();
    assertEquals(builder, builder.outputParams(null));
    assertEquals(builder, builder.filterParams(SamFilterParams.builder().create()));
    assertEquals(builder, builder.readGroupLabels(null));
    assertEquals(builder, builder.readGroupStatistics(null));
  }

}
