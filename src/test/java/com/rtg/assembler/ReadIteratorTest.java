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

import java.io.IOException;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import com.rtg.mode.DnaUtils;
import com.rtg.reader.ReaderTestUtils;
import com.rtg.reader.SequencesReader;
import com.rtg.util.diagnostic.Diagnostic;

import junit.framework.TestCase;

/**
 */
public class ReadIteratorTest extends TestCase {
  public void test() throws IOException {
    final String[] expected = {
        "ACGGT"
        , "TTG"
        , "TATG"
    };
    check(expected);
  }
  public void testEmpty() throws IOException {
    final String[] expected = {};
    check(expected);
  }

  public void testSingle() throws IOException {
    final String[] expected = {"A"};
    check(expected);
  }

  static void listEquals(String[] a, List<byte[]> b) {
    assertEquals("expected <" + Arrays.toString(a) + "> but was <" + ReadPairSourceTest.fragementToString(b) + ">", a.length, b.size());
    for (int i = 0; i < a.length; ++i) {
      assertEquals(a[i], DnaUtils.bytesToSequenceIncCG(b.get(i)));
    }
  }
  private void check(String[] expected) throws IOException {
    final String fasta = ReaderTestUtils.fasta(expected);
    final SequencesReader reader = ReaderTestUtils.getReaderDnaMemory(fasta);
    final AsyncReadPool pool = new AsyncReadPool("ReadIteratorTest", Collections.singletonList(new ReadPairSource(reader)));
    final ReadIterator iterator = new ReadIterator(pool.sources().get(0));

    for (final String s : expected) {
      assertTrue(iterator.hasNext());

      listEquals(new String[] {s}, iterator.next());
    }
    if (iterator.hasNext()) {
      fail("Additional sequence at end: " + ReadPairSourceTest.fragementToString(iterator.next()));
    }
    pool.terminate();
  }
  public void testPaired() throws IOException {
    Diagnostic.setLogStream();
    final String[] left = {
        "ACGGGT"
        , "ACCGGT"
        , "AAAAAA"
    };
    final String[] right = {
        "AAGAGAAGA"
        , "CCCCCC"
        , "CGCGAAGCGC"
    };
    final String leftFasta = ReaderTestUtils.fasta(left);
    final String rightFasta = ReaderTestUtils.fasta(right);
    final ReadPairSource reader = new ReadPairSource(ReaderTestUtils.getReaderDnaMemory(leftFasta), ReaderTestUtils.getReaderDnaMemory(rightFasta));

    for (int j = 0; j < 2; ++j) {
      reader.reset();
      final AsyncReadPool pool = new AsyncReadPool("ReadIteratorTest", Collections.singletonList(reader));
      final ReadIterator iterator = new ReadIterator(pool.sources().get(0));
      int i = 0;
      while (iterator.hasNext()) {
        listEquals(new String[] {left[i], right[i]}, iterator.next());
        ++i;
      }
      assertEquals(3, i);
      pool.terminate();
    }
  }
}
