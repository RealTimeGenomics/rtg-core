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
package com.rtg.index.hash.ngs;


import com.rtg.util.TestUtils;

import junit.framework.TestCase;
import junit.framework.TestSuite;

/**
 */
public class ReadEncoderTest extends TestCase {

  public final void test() {
    TestUtils.testPseudoEnum(ReadEncoder.class, "[PAIRED_FIRST, PAIRED_SECOND, SINGLE_END]");
  }

  public final void testFirst() {
    assertEquals(0, ReadEncoder.PAIRED_FIRST.encode(0));
    assertEquals(1, ReadEncoder.PAIRED_SECOND.encode(0));

    assertEquals(6, ReadEncoder.PAIRED_FIRST.encode(3));
    assertEquals(7, ReadEncoder.PAIRED_SECOND.encode(3));

    assertEquals(0, ReadEncoder.SINGLE_END.encode(0));
    assertEquals(3, ReadEncoder.SINGLE_END.encode(3));
  }

  public static void main(final String[] args) {
    junit.textui.TestRunner.run(new TestSuite(ReadEncoderTest.class));
  }

}

