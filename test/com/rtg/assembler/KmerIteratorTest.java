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

import junit.framework.TestCase;

/**
 */
public class KmerIteratorTest extends TestCase {
  public void test() {
    KmerIterator iterator = new KmerIterator(new byte[] {1, 2, 3, 1, 0, 2, 3, 4, 1}, StringKmer.factory(), 4);
    String[] expected = {
        "ACGA"
        , "CGTA"
    };
    check(iterator, expected);
  }

  public void testEmpty() {
    KmerIterator iterator = new KmerIterator(new byte[] {}, StringKmer.factory(), 4);
    String[] expected = {
    };
    check(iterator, expected);
  }
  public void testShort() {
    KmerIterator iterator = new KmerIterator(new byte[] {1, 2, 3}, StringKmer.factory(), 4);
    String[] expected = {
    };
    check(iterator, expected);
  }
  public void testToManyNs() {
    KmerIterator iterator = new KmerIterator(new byte[] {1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0}, StringKmer.factory(), 4);
    String[] expected = {
    };
    check(iterator, expected);
  }
  public void testNStart() {
    KmerIterator iterator = new KmerIterator(new byte[] {0, 1, 2, 3, 4, 2, 2, 3, 2, 1, 2, 3, 2}, StringKmer.factory(), 4);
    String[] expected = {
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
