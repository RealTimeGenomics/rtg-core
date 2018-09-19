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

import com.rtg.assembler.graph.Contig;
import com.rtg.mode.DNARange;

import junit.framework.TestCase;

/**
 */
public class ContigStringTest extends TestCase {

  public void test() {
    final Contig contig = new ContigString("ACGTN");
    assertEquals(5, contig.length());
    assertEquals(DNARange.A, contig.nt(0));
    assertEquals(DNARange.C, contig.nt(1));
    assertEquals(DNARange.G, contig.nt(2));
    assertEquals(DNARange.T, contig.nt(3));
    assertEquals(DNARange.N, contig.nt(4));

    try {
      contig.nt(-1);
      fail();
    } catch (final RuntimeException e) {
      // expected
    }
    try {
      contig.nt(5);
      fail();
    } catch (final RuntimeException e) {
      // expected
    }
  }

  public void test0() {
    final Contig contig = new ContigString("");
    assertEquals(0, contig.length());
    try {
      contig.nt(-1);
      fail();
    } catch (final RuntimeException e) {
      // expected
    }
    try {
      contig.nt(0);
      fail();
    } catch (final RuntimeException e) {
      // expected
    }
  }

  public void testContigSequenceString() {
    final String nt = "ACGTN";
    final Contig contig = new ContigString(nt);
    assertEquals(nt, ContigString.contigSequenceString(contig));
  }
}
