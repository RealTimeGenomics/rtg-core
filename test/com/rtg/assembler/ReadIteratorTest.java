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
    for (int i = 0; i < a.length; i++) {
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

    for (int j = 0; j < 2; j++) {
      reader.reset();
      final AsyncReadPool pool = new AsyncReadPool("ReadIteratorTest", Collections.singletonList(reader));
      final ReadIterator iterator = new ReadIterator(pool.sources().get(0));
      int i = 0;
      while (iterator.hasNext()) {
        listEquals(new String[] {left[i], right[i]}, iterator.next());
        i++;
      }
      assertEquals(3, i);
      pool.terminate();
    }
  }
}
