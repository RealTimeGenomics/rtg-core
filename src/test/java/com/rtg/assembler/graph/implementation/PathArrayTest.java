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
