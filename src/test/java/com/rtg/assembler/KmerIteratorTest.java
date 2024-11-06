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

package com.rtg.assembler;

import junit.framework.TestCase;

/**
 */
public class KmerIteratorTest extends TestCase {
  public void test() {
    final KmerIterator iterator = new KmerIterator(new byte[] {1, 2, 3, 1, 0, 2, 3, 4, 1}, StringKmer.factory(), 4);
    final String[] expected = {
        "ACGA"
        , "CGTA"
    };
    check(iterator, expected);
  }

  public void testEmpty() {
    final KmerIterator iterator = new KmerIterator(new byte[] {}, StringKmer.factory(), 4);
    final String[] expected = {
    };
    check(iterator, expected);
  }
  public void testShort() {
    final KmerIterator iterator = new KmerIterator(new byte[] {1, 2, 3}, StringKmer.factory(), 4);
    final String[] expected = {
    };
    check(iterator, expected);
  }
  public void testToManyNs() {
    final KmerIterator iterator = new KmerIterator(new byte[] {1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0}, StringKmer.factory(), 4);
    final String[] expected = {
    };
    check(iterator, expected);
  }
  public void testNStart() {
    final KmerIterator iterator = new KmerIterator(new byte[] {0, 1, 2, 3, 4, 2, 2, 3, 2, 1, 2, 3, 2}, StringKmer.factory(), 4);
    final String[] expected = {
        "ACGT"
        , "CGTC"
        , "GTCC"
        , "TCCG"
        , "CCGC"
        , "CGCA"
        , "GCAC"
        , "CACG"
        , "ACGC"
    };
    check(iterator, expected);
  }

  private void check(KmerIterator iterator, String[] expected) {
    for (String anExpected : expected) {
      assertTrue(anExpected, iterator.hasNext());
      assertEquals(anExpected, iterator.next().toString());
    }
    assertFalse(iterator.hasNext());
  }
}
