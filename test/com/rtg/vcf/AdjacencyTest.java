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

package com.rtg.vcf;

import java.io.IOException;

import junit.framework.TestCase;

/**
 */
public class AdjacencyTest extends TestCase {

  public void test() {
    final Adjacency adj = new Adjacency("chr1", 10, true, "ACGT", "chr2", 22, false);
    assertEquals("chr1", adj.thisChromosome());
    assertEquals(10,     adj.thisEndPos());
    assertEquals("ACGT", adj.thisBases());
    assertEquals(true,   adj.thisIsForward());

    assertEquals("chr2", adj.mateChromosome());
    assertEquals(22,     adj.mateStartPos());
    assertEquals(false,  adj.mateIsForward());
  }

  public void testParseFwdFwd() throws IOException {
    final Adjacency adj = Adjacency.parseAdjacency("chr1", 12, "A[chr2:9[");
    assertEquals("chr1", adj.thisChromosome());
    assertEquals(12,     adj.thisEndPos());
    assertEquals("A",    adj.thisBases());
    assertEquals(true,   adj.thisIsForward());

    assertEquals("chr2", adj.mateChromosome());
    assertEquals(9,     adj.mateStartPos());
    assertEquals(true,  adj.mateIsForward());
  }

  public void testParseFwdRev() throws IOException {
    final Adjacency adj = Adjacency.parseAdjacency("chr1", 12, "A]chr2:9]");
    assertEquals("chr1", adj.thisChromosome());
    assertEquals(12,     adj.thisEndPos());
    assertEquals("A",    adj.thisBases());
    assertEquals(true,   adj.thisIsForward());

    assertEquals("chr2", adj.mateChromosome());
    assertEquals(9,     adj.mateStartPos());
    assertEquals(false,  adj.mateIsForward());
  }

  public void testParseRevFwd() throws IOException {
    final Adjacency adj = Adjacency.parseAdjacency("chr1", 12, "[chr2:9[CG");
    assertEquals("chr1", adj.thisChromosome());
    assertEquals(12,     adj.thisEndPos());
    assertEquals("CG",   adj.thisBases());
    assertEquals(false,   adj.thisIsForward());

    assertEquals("chr2", adj.mateChromosome());
    assertEquals(9,     adj.mateStartPos());
    assertEquals(true,  adj.mateIsForward());
  }

  public void testParseRevRev() throws IOException {
    final Adjacency adj = Adjacency.parseAdjacency("chr1", 12, "]chr2:9]CG");
    assertEquals("chr1", adj.thisChromosome());
    assertEquals(12,     adj.thisEndPos());
    assertEquals("CG",   adj.thisBases());
    assertEquals(false,   adj.thisIsForward());

    assertEquals("chr2", adj.mateChromosome());
    assertEquals(9,     adj.mateStartPos());
    assertEquals(false,  adj.mateIsForward());
  }

  public void testBadParse() throws IOException {
    assertNull(Adjacency.parseAdjacency("chr1", 12, "C"));
    assertNull(Adjacency.parseAdjacency("chr1", 12, "<DEL>"));
    assertNull(Adjacency.parseAdjacency("chr1", 12, "<DUP:TANDEM>"));
  }

  public void testExceptionParse() {
    checkIoException("C[]");
    checkIoException("C[chr:32");
    checkIoException("C(chr:32[");
    checkIoException("C[chr:32[C");
    checkIoException("(chr:32[C");
    checkIoException("[chr:32(C");
  }

  protected void checkIoException(final String call) {
    try {
      Adjacency.parseAdjacency("chr1", 12, call);
      fail("Excepted exception when parsing: ");
    } catch (IOException ex) {
      assertEquals("Invalid .vcf adjacency call: " + call, ex.getMessage());
    }
  }

  public void testCompareTo() {
    // in strictly sorted order
    final Adjacency[] a = {
        new Adjacency("chr1", 10, true),
    new Adjacency("chr1", 10, false),
    new Adjacency("chr1", 11, true),
    new Adjacency("chr2", 10, true),
    };
    // check all pairs
    for (int i = 0; i < a.length; i++) {
      for (int j = 0; j < a.length; j++) {
        final int expect = Integer.signum(i - j);
        assertEquals(expect, a[i].compareTo(a[j]));
        if (expect == 0) {
          assertEquals(a[i].hashCode(), a[j].hashCode());
        }
      }
    }
  }
}
