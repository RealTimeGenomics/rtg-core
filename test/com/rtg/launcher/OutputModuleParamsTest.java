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

package com.rtg.launcher;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;

import com.rtg.launcher.OutputModuleParams.OutputModuleParamsBuilder;
import com.rtg.launcher.OutputModuleParamsTest.MockOutputModuleParams.MockOutputModuleParamsBuilder;
import com.rtg.util.StringUtils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.FileUtils;
import com.rtg.util.test.FileHelper;
import com.rtg.util.test.params.TestParams;

import junit.framework.TestCase;

/**
 */
public class OutputModuleParamsTest extends TestCase {

  static final class MockOutputModuleParams extends OutputModuleParams {

    static MockOutputModuleParamsBuilder builder() {
      return new MockOutputModuleParamsBuilder();
    }

    static final class MockOutputModuleParamsBuilder extends OutputModuleParamsBuilder<MockOutputModuleParamsBuilder> {
      @Override
      protected MockOutputModuleParamsBuilder self() {
        return this;
      }

      public MockOutputModuleParams create() {
        return new MockOutputModuleParams(this);
      }
    }

    public MockOutputModuleParams(MockOutputModuleParamsBuilder builder) {
      super(builder);
    }
  }

  public void testOmnes() {
    new TestParams(OutputModuleParams.class, OutputModuleParamsBuilder.class).check();
  }

  public void testDefaultOutputModuleParams() throws IOException {
    Diagnostic.setLogStream();
    final File tempDir = FileUtils.createTempDir("outputmoduleparams", "test");
    try {
      final MockOutputModuleParams def = MockOutputModuleParams.builder().outputParams(new OutputParams(tempDir, false, true)).create();
      assertNotNull(def.outputParams());
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
      assertEquals("OutputParams output directory=" + tempDir.getPath() + " progress=" + false + " zip=" + true + StringUtils.LS, def.toString());
    } finally {
      assertTrue(FileHelper.deleteAll(tempDir));
    }
  }

  public void testOutputModuleParamsBuilder() {
    final MockOutputModuleParamsBuilder builder = MockOutputModuleParams.builder();
    assertEquals(builder, builder.outputParams(null));
  }

}
