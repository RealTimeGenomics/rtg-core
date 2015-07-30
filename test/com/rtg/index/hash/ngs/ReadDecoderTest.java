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

/**
 */
public class ReadDecoderTest extends TestCase {

  public final void test() {
    TestUtils.testPseudoEnum(ReadDecoder.class, "[PAIRED_END, SINGLE_END]");
  }

  public final void testFirst() {

    assertEquals(0, ReadDecoder.PAIRED_END.decode(0));
    assertEquals(3, ReadDecoder.PAIRED_END.decode(7));
    assertEquals(2, ReadDecoder.PAIRED_END.decode(4));
    assertTrue(ReadDecoder.PAIRED_END.isFirst(0));
    assertTrue(ReadDecoder.PAIRED_END.isFirst(4));
    assertFalse(ReadDecoder.PAIRED_END.isFirst(7));

    assertEquals(0, ReadDecoder.SINGLE_END.decode(0));
    assertEquals(7, ReadDecoder.SINGLE_END.decode(7));
    assertTrue(ReadDecoder.SINGLE_END.isFirst(0));
    assertTrue(ReadDecoder.SINGLE_END.isFirst(7));
  }
}

