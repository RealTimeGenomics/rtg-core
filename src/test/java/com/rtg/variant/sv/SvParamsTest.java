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
