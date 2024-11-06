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
