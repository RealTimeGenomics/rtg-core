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

package com.rtg.sam;

import java.io.File;
import java.io.IOException;

import com.rtg.launcher.MockReaderParams;
import com.rtg.launcher.MockSequenceParams;
import com.rtg.launcher.OutputParams;
import com.rtg.mode.SequenceMode;
import com.rtg.sam.MappedParams.MappedParamsBuilder;
import com.rtg.sam.MappedParamsTest.MockMappedParams.MockMappedParamsBuilder;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.FileUtils;
import com.rtg.util.test.FileHelper;
import com.rtg.util.test.params.TestParams;

import junit.framework.TestCase;

/**
 */
public class MappedParamsTest extends TestCase {

  static final class MockMappedParams extends MappedParams {

    static MockMappedParamsBuilder builder() {
      return new MockMappedParamsBuilder();
    }

    static final class MockMappedParamsBuilder extends MappedParamsBuilder<MockMappedParamsBuilder> {
      @Override
      protected MockMappedParamsBuilder self() {
        return this;
      }

      public MockMappedParams create() {
        return new MockMappedParams(this);
      }
    }

    public MockMappedParams(MockMappedParamsBuilder builder) {
      super(builder);
    }
  }

  public void testOmnes() {
    new TestParams(MappedParams.class, MappedParamsBuilder.class).check();
  }

  public void testDefaultMappedParams() throws IOException {
    Diagnostic.setLogStream();
    final File tempDir = FileUtils.createTempDir("mappedparams", "test");
    try {
      MockMappedParams def = MockMappedParams.builder().outputParams(new OutputParams(tempDir, false, true)).create();
      assertNotNull(def.filterParams());
      assertNull(def.genome());
      def.close();
      def = MockMappedParams.builder().outputParams(new OutputParams(tempDir, false, true)).genome(new MockSequenceParams(new MockReaderParams(1, 1, SequenceMode.BIDIRECTIONAL))).create();
      assertNotNull(def.genome());
      def.close();
    } finally {
      assertTrue(FileHelper.deleteAll(tempDir));
    }
  }

  public void testMappedParamsBuilder() {
    final MockMappedParamsBuilder builder = MockMappedParams.builder();
    assertEquals(builder, builder.genome(null));
    assertEquals(builder, builder.filterParams(null));
  }
}
