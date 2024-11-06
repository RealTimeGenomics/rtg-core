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
package com.rtg.index.similarity;

import com.rtg.util.StringUtils;

import junit.framework.TestCase;

/**
 */
public class SimilaritySorterTest extends TestCase {

  public void testSmall() {
    final SimilaritySorter sorter = new SimilaritySorter(2, false);
    sorter.globalIntegrity();

    sorter.add(0);
    sorter.add(0);
    sorter.globalIntegrity();
    final String expUnsorted = ""
      + "SimilaritySorter 2:2 unsorted" + StringUtils.LS
      + "0 0 " + StringUtils.LS
      ;
    assertEquals(expUnsorted, sorter.toString());

    final SimilarityMatrix matrix = new SimilarityMatrix(1);
    sorter.similarity(matrix);
    sorter.globalIntegrity();
    final String expSorted = ""
      + "SimilaritySorter 2:2 sorted" + StringUtils.LS
      + "0 0 " + StringUtils.LS
      + "duplicates removed 1" + StringUtils.LS
      + "0 " + StringUtils.LS
      + "2 " + StringUtils.LS
      ;
    assertEquals(expSorted, sorter.toString());
    //System.err.println(matrix.toString());
    assertEquals(4, matrix.get(0, 0), 1E-8);

    sorter.reset();
    sorter.globalIntegrity();
  }

  /** Covers a jumble case */
  public void test1() {
    final SimilaritySorter sorter = new SimilaritySorter(2, false);
    sorter.globalIntegrity();

    sorter.add(1);
    sorter.add(1);
    sorter.globalIntegrity();
    final String expUnsorted = ""
      + "SimilaritySorter 2:2 unsorted" + StringUtils.LS
      + "1 1 " + StringUtils.LS
      ;
    assertEquals(expUnsorted, sorter.toString());

    final SimilarityMatrix matrix = new SimilarityMatrix(2);
    sorter.similarity(matrix);
    sorter.globalIntegrity();
    final String expSorted = ""
      + "SimilaritySorter 2:2 sorted" + StringUtils.LS
      + "1 1 " + StringUtils.LS
      + "duplicates removed 1" + StringUtils.LS
      + "1 " + StringUtils.LS
      + "2 " + StringUtils.LS
      ;
    assertEquals(expSorted, sorter.toString());
    //System.err.println(matrix.toString());
    assertEquals(0, matrix.get(0, 0), 1E-8);
    assertEquals(0, matrix.get(0, 1), 1E-8);
    assertEquals(4, matrix.get(1, 1), 1E-8);

    sorter.reset();
    sorter.globalIntegrity();
  }

  public void testReset() {
    final SimilaritySorter sorter = new SimilaritySorter(10, false);
    sorter.globalIntegrity();

    sorter.add(1);
    sorter.add(2);
    sorter.add(3);
    sorter.add(0);
    sorter.add(1);
    sorter.add(2);
    sorter.globalIntegrity();
    final String expUnsorted = ""
      + "SimilaritySorter 6:10 unsorted" + StringUtils.LS
      + "1 2 3 0 1 2 " + StringUtils.LS
      ;
    assertEquals(expUnsorted, sorter.toString());

    final SimilarityMatrix matrix = new SimilarityMatrix(4);
    sorter.similarity(matrix);
    sorter.globalIntegrity();
    final String expSorted = ""
      + "SimilaritySorter 6:10 sorted" + StringUtils.LS
      + "0 1 1 2 2 3 " + StringUtils.LS
      + "duplicates removed 4" + StringUtils.LS
      + "0 1 2 3 " + StringUtils.LS
      + "1 2 2 1 " + StringUtils.LS
      ;
    assertEquals(expSorted, sorter.toString());
    //System.err.println(matrix.toString());
    assertEquals(1, matrix.get(0, 0), 1E-8);
    assertEquals(2, matrix.get(0, 1), 1E-8);
    assertEquals(2, matrix.get(0, 2), 1E-8);
    assertEquals(1, matrix.get(0, 3), 1E-8);

    assertEquals(4, matrix.get(1, 1), 1E-8);
    assertEquals(4, matrix.get(1, 2), 1E-8);
    assertEquals(2, matrix.get(1, 3), 1E-8);

    assertEquals(4, matrix.get(2, 2), 1E-8);
    assertEquals(2, matrix.get(2, 3), 1E-8);

    assertEquals(1, matrix.get(3, 3), 1E-8);

    sorter.reset();
    sorter.globalIntegrity();

    sorter.add(1);
    sorter.add(2);

    sorter.globalIntegrity();
    final String expUnsorted2 = ""
      + "SimilaritySorter 2:10 unsorted" + StringUtils.LS
      + "1 2 " + StringUtils.LS
      ;
    assertEquals(expUnsorted2, sorter.toString());

    final SimilarityMatrix matrix2 = new SimilarityMatrix(4);
    sorter.similarity(matrix2);
    sorter.globalIntegrity();
    final String expSorted2 = ""
      + "SimilaritySorter 2:10 sorted" + StringUtils.LS
      + "1 2 " + StringUtils.LS
      + "duplicates removed 2" + StringUtils.LS
      + "1 2 " + StringUtils.LS
      + "1 1 " + StringUtils.LS
      ;
    assertEquals(expSorted2, sorter.toString());

  }

}
