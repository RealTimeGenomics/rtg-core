/*
 * Copyright (c) 2017. Real Time Genomics Limited.
 *
 * Use of this source code is bound by the Real Time Genomics Limited Software Licence Agreement
 * for Academic Non-commercial Research Purposes only.
 *
 * If you did not receive a license accompanying this file, a copy must first be obtained by email
 * from support@realtimegenomics.com.  On downloading, using and/or continuing to use this source
 * code you accept the terms of that license agreement and any amendments to those terms that may
 * be made from time to time by Real Time Genomics Limited.
 */

package com.rtg.reader;

import static org.junit.Assert.assertEquals;

import java.io.StringWriter;

import org.junit.Test;

/**
 */
public class AsyncFastqPairWriterTest {
  @Test
  public void testWrite() {
    final StringWriter leftOut = new StringWriter();
    final StringWriter rightOut = new StringWriter();
    final FastqWriter left = new FastqWriter(leftOut);
    final FastqWriter right = new FastqWriter(rightOut);
    try (final AsyncFastqPairWriter writer = new AsyncFastqPairWriter(left, right)) {
      writer.write(new FastqPair(FastqSequenceTest.getFastq("foo", "ACGT"), FastqSequenceTest.getFastq("foo", "TTTT")));
    }
    assertEquals("@foo\nACGT\n+\naaaa\n", leftOut.toString());
    assertEquals("@foo\nTTTT\n+\naaaa\n", rightOut.toString());
  }

}