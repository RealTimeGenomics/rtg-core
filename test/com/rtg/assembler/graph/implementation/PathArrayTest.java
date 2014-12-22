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

package com.rtg.assembler.graph.implementation;

import java.util.Arrays;

import com.rtg.assembler.graph.Path;

import junit.framework.TestCase;

/**
 */
public class PathArrayTest extends TestCase {

  public void test() {
    final PathArray path = new PathArray(+1, -2, +3);
    path.globalIntegrity();
    assertEquals(3, path.length());
    assertEquals(+1, path.contig(0));
    assertEquals(-2, path.contig(1));
    assertEquals(+3, path.contig(2));

    assertEquals(0, path.index(+1));
    assertEquals(1, path.index(-2));
    assertEquals(2, path.index(+3));

    assertEquals(-1, path.index(-1));
    assertEquals(-1, path.index(+4));

    try {
      path.index(0);
      fail();
    } catch (final IllegalArgumentException e) {
      // expected
    }
  }

  public void test1() {
    final PathArray path = new PathArray(+1);
    try {
      path.globalIntegrity();
      fail();
    } catch (final Exception e) {
      // expected
    }
  }

  public void test2() {
    final PathArray path = new PathArray(+1, -2);
    path.globalIntegrity();
    assertEquals(2, path.length());
    assertEquals(+1, path.contig(0));
    assertEquals(-2, path.contig(1));
    try {
      path.contig(-1);
      fail();
    } catch (final RuntimeException e) {
      // TODO: handle exception
    }
    try {
      path.contig(2);
      fail();
    } catch (final RuntimeException e) {
      // TODO: handle exception
    }

    assertEquals(0, path.index(+1));
    assertEquals(1, path.index(-2));

    assertEquals(-1, path.index(-1));
    assertEquals(-1, path.index(+4));

    try {
      path.index(0);
      fail();
    } catch (final IllegalArgumentException e) {
      // expected
    }
  }
  public void testFromList() {
    final Path path = new PathArray(Arrays.asList(1L, 5L, 9L));
    assertEquals(3, path.length());
    assertEquals(1L, path.contig(0));
    assertEquals(5L, path.contig(1));
    assertEquals(9L, path.contig(2));
  }


}
