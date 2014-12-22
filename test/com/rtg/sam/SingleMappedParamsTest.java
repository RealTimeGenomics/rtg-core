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
import java.util.ArrayList;

import com.rtg.sam.SingleMappedParams.SingleMappedParamsBuilder;
import com.rtg.util.test.params.TestParams;

import junit.framework.TestCase;

/**
 */
public class SingleMappedParamsTest extends TestCase {

  /**
   * Mock SingleMappedParamsBuilder class
   */
  public static final class MockSingleMappedParamsBuilder extends SingleMappedParamsBuilder<MockSingleMappedParamsBuilder> {
    @Override
    protected MockSingleMappedParamsBuilder self() {
      return this;
    }
  }

  static final class MockSingleMappedParams extends SingleMappedParams {

    MockSingleMappedParams(MockSingleMappedParamsBuilder builder) {
      super(builder);
    }
  }

  public void testSingleMappedParamsBuilder() {
    final MockSingleMappedParams dummy = new MockSingleMappedParams(new MockSingleMappedParamsBuilder().mapped(new ArrayList<File>()).ioThreads(2).execThreads(3));
    assertEquals(2, dummy.ioThreads());
    assertEquals(3, dummy.execThreads());
    assertEquals(0, dummy.mapped().size());
  }

  public void testSingleMappedParamsDefaults() {
    final MockSingleMappedParams dummy = new MockSingleMappedParams(new MockSingleMappedParamsBuilder());
    assertEquals(1, dummy.ioThreads());
    assertEquals(1, dummy.execThreads());
    assertNull(dummy.mapped());
  }

  public void testOmnes() {
    new TestParams(SingleMappedParams.class, SingleMappedParamsBuilder.class).check();
  }
}
