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
package com.rtg.ngs;


import static com.rtg.util.StringUtils.LS;

import java.io.ByteArrayOutputStream;
import java.io.IOException;

import com.rtg.launcher.HashingRegion;

import junit.framework.TestCase;

/**
 * Test default output processor
 *
 *
 *
 */
public class DefaultOutputProcessorTest extends TestCase {

  public void testProcessor() throws IOException {
    final ByteArrayOutputStream b = new ByteArrayOutputStream();
    try {
      final DefaultOutputProcessor dop = new DefaultOutputProcessor(NgsParams.builder().outputParams(new NgsTestUtils.OverriddenNgsOutputParams(NgsTestUtils.OverriddenNgsOutputParams.builder().outStream(b))).create());
      dop.process(1234, "R", 123456, 30000, 2, 0);
      dop.finish();
    } finally {
      b.close();
    }
    assertEquals("#template-id\tframe\tread-id\ttemplate-start\tscore\tscore-indel" + LS + "1234\tR\t123456\t30001\t2\t0" + LS, b.toString());
  }

  public void testThreadClone() throws IOException {
    final ByteArrayOutputStream b = new ByteArrayOutputStream();
    final DefaultOutputProcessor dop = new DefaultOutputProcessor(NgsParams.builder().outputParams(new NgsTestUtils.OverriddenNgsOutputParams(NgsTestUtils.OverriddenNgsOutputParams.builder().outStream(b))).create());
    try {
      dop.threadClone(HashingRegion.NONE);
      fail();
    } catch (final UnsupportedOperationException e) {
      assertEquals("DefaultOutputProcessor is not thread-safe", e.getMessage());
    }
  }
}

