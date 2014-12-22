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
