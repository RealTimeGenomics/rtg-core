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

package com.rtg.reader;

import static org.junit.Assert.assertEquals;

import java.io.StringWriter;
import java.util.Arrays;

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

  @Test
  public void testAccept() {
    final StringWriter leftOut = new StringWriter();
    final StringWriter rightOut = new StringWriter();
    final FastqWriter left = new FastqWriter(leftOut);
    final FastqWriter right = new FastqWriter(rightOut);
    try (final AsyncFastqPairWriter writer = new AsyncFastqPairWriter(left, right)) {
      writer.accept(Arrays.asList(
        new FastqPair(FastqSequenceTest.getFastq("foo", "ACGT"), FastqSequenceTest.getFastq("foo", "TTTT")),
        new FastqPair(FastqSequenceTest.getFastq("bar", "CCCC"), FastqSequenceTest.getFastq("bar", "TTAA"))
      ));
    }
    assertEquals("@foo\nACGT\n+\naaaa\n@bar\nCCCC\n+\naaaa\n", leftOut.toString());
    assertEquals("@foo\nTTTT\n+\naaaa\n@bar\nTTAA\n+\naaaa\n", rightOut.toString());
  }
}
