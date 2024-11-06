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

package com.rtg.assembler;

import com.rtg.assembler.graph.implementation.ContigString;
import com.rtg.util.diagnostic.Diagnostic;

/**
 */
public class ByteKmerTest extends AbstractKmerTest {

  @Override
  public void setUp() {
    Diagnostic.setLogStream();
  }

  @Override
  protected Kmer kmer(String str) {
    return new ByteKmer(str);
  }

  public void testMinimalKmer() {
    ByteKmer kmer = new ByteKmer("ACGT");
    assertEquals("ACGT", kmer.minimalKmer().toString());
    kmer = new ByteKmer("TTTT");
    assertEquals("AAAA", kmer.minimalKmer().toString());
    kmer = new ByteKmer("TTTG");
    assertEquals("CAAA", kmer.minimalKmer().toString());

    kmer = new ByteKmer("TTTGA");
    assertEquals("TCAAA", kmer.minimalKmer().toString());
    kmer = new ByteKmer("CGGC");
    assertEquals("CGGC", kmer.minimalKmer().toString());
    kmer = new ByteKmer("GCCG");
    assertEquals("CGGC", kmer.minimalKmer().toString());
  }

  public void testSetGetPos() {
    final byte[] bytes = new byte[5];
    assertEquals(0, ByteKmer.getPos(bytes, 19));
    ByteKmer.setPos(bytes, 19, 3);
    assertEquals(3, ByteKmer.getPos(bytes, 19));
    ByteKmer.setPos(bytes, 19, 2);
    assertEquals(2, ByteKmer.getPos(bytes, 19));

    for (int i = 0; i < 20; ++i) {
      ByteKmer.setPos(bytes, i, 3);
    }
    ByteKmer.setPos(bytes, 19, 0);
    assertEquals(0, ByteKmer.getPos(bytes, 19));
    ByteKmer.setPos(bytes, 5, 2);
    assertEquals(2, ByteKmer.getPos(bytes, 5));

    ByteKmer.setPos(bytes, 4, 2);
    assertEquals(2, ByteKmer.getPos(bytes, 4));
  }

  public void testByteFactory() {
    final KmerFactory fact = ByteKmer.factory();
    assertEquals("ACGT", fact.make(new byte[] {1, 1, 1, 1, 2, 3, 4, 2, 2, 2, 2},  3,  7).toString());
    assertEquals("AGAAC", fact.make(new byte[] {1, 3, 1, 1, 2, 3, 4, 2, 2, 2, 2},  0,  5).toString());
    assertEquals("TCCGC", fact.make(new byte[] {1, 3, 1, 1, 2, 3, 4, 2, 2, 3, 2},  6,  11).toString());

  }
  public void testContigFactory() {
    final KmerFactory fact = ByteKmer.factory();
    assertEquals("ACGT", fact.make(new ContigString("AAAACGTCCCC"), 3, 7).toString());
    assertEquals("AGAAC", fact.make(new ContigString("AGAACTGTAGT"),  0,  5).toString());
    assertEquals("TCCGC", fact.make(new ContigString("AGAACGTCCGC"),  6,  11).toString());
  }

  public void testGetPos() {
    final byte[] bs = {-1, 15, 0};
    assertEquals(3, ByteKmer.getPos(bs, 0));
    assertEquals(3, ByteKmer.getPos(bs, 3));
    assertEquals(3, ByteKmer.getPos(bs, 4));
    assertEquals(0, ByteKmer.getPos(bs, 7));
    assertEquals(0, ByteKmer.getPos(bs, 8));
    assertEquals(0, ByteKmer.getPos(bs, 11));

    try {
      ByteKmer.getPos(bs, -1);
      fail();
    } catch (final Exception e) {
      assertEquals("pos=-1", e.getMessage());
    }
    try {
      ByteKmer.getPos(bs, 12);
      fail();
    } catch (final Exception e) {
      assertEquals("pos=12", e.getMessage());
    }
  }

  public void testSetPos() {
    final byte[] bs = {-1, 15, 0};
    try {
      ByteKmer.setPos(bs, -1, 0);
      fail();
    } catch (final Exception e) {
      assertEquals("pos=-1", e.getMessage());
    }
    try {
      ByteKmer.setPos(bs, 12, 0);
      fail();
    } catch (final Exception e) {
      assertEquals("pos=12", e.getMessage());
    }
  }
}
