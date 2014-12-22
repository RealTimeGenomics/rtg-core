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
package com.rtg.reader;

import junit.framework.TestCase;

/**
 */
public class PointerFileLookupTest extends TestCase {
  public void testSpanning() {
    int[][] pointers = {
            new int[] {0, 35, 60}, //sequences 0 and start of 1
            new int[] {20, 19},    //end of 1, sequence 2
            new int[] {40},        //cont
            new int[] {3}          //fin
    };
    PointerFileLookup p = PointerFileLookup.generateLookup(pointers);
    assertEquals(0, p.lookup(0));
    assertEquals(0, p.lookup(1));
    assertEquals(1, p.lookup(2));
    assertEquals(3, p.lookup(3));
    assertEquals(0, p.startSeq(0));
    assertEquals(2, p.startSeq(1));
    assertEquals(3, p.startSeq(2));
    assertEquals(3, p.startSeq(3));
  }

  public void testNormal() {
    int[][] pointers = {
            new int[] {0, 35, 60, 90},
            new int[] {20, 19, 40, 60},
            new int[] {40, 50, 60, 70},
            new int[] {30, 40, 100, 200}
    };
    PointerFileLookup p = PointerFileLookup.generateLookup(pointers);
    assertEquals(0, p.lookup(0));
    assertEquals(0, p.lookup(1));
    assertEquals(0, p.lookup(2));
    assertEquals(1, p.lookup(3));
    assertEquals(1, p.lookup(4));
    assertEquals(1, p.lookup(5));
    assertEquals(2, p.lookup(6));
    assertEquals(2, p.lookup(7));
    assertEquals(2, p.lookup(8));
    assertEquals(3, p.lookup(9));
    assertEquals(3, p.lookup(10));
    assertEquals(3, p.lookup(11));
    assertEquals(3, p.lookup(12));
    assertEquals(0, p.startSeq(0));
    assertEquals(3, p.startSeq(1));
    assertEquals(6, p.startSeq(2));
    assertEquals(9, p.startSeq(3));
  }

  public void testSingleFile() {
    int[][] pointers = {
            new int[] {0, 35, 60, 90, 100}, //sequences 0, 1, 2, 3
    };
    PointerFileLookup p = PointerFileLookup.generateLookup(pointers);
    assertEquals(0, p.lookup(0));
    assertEquals(0, p.lookup(1));
    assertEquals(0, p.lookup(2));
    assertEquals(0, p.lookup(3));
    assertEquals(0, p.lookup(4));
    assertEquals(0, p.startSeq(0));
  }
  public void testDoubleFile() {
    int[][] pointers = {
            new int[] {0, 35, 60, 90, 100}, //sequences 0, 1, 2, 3
            new int[] {0, 35, 60, 90, 100}, //sequences 4, 5, 6, 7
    };
    PointerFileLookup p = PointerFileLookup.generateLookup(pointers);
    assertEquals(0, p.lookup(0));
    assertEquals(0, p.lookup(1));
    assertEquals(0, p.lookup(2));
    assertEquals(0, p.lookup(3));
    assertEquals(1, p.lookup(4));
    assertEquals(1, p.lookup(5));
    assertEquals(1, p.lookup(6));
    assertEquals(1, p.lookup(7));
    assertEquals(1, p.lookup(8));
    assertEquals(0, p.startSeq(0));
    assertEquals(4, p.startSeq(1));
  }
}
