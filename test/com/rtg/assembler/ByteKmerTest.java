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
      assertEquals("-1", e.getMessage());
    }
    try {
      ByteKmer.getPos(bs, 12);
      fail();
    } catch (final Exception e) {
      assertEquals("12", e.getMessage());
    }
  }

  public void testSetPos() {
    final byte[] bs = {-1, 15, 0};
    try {
      ByteKmer.setPos(bs, -1, 0);
      fail();
    } catch (final Exception e) {
      assertEquals("-1", e.getMessage());
    }
    try {
      ByteKmer.setPos(bs, 12, 0);
      fail();
    } catch (final Exception e) {
      assertEquals("12", e.getMessage());
    }
  }
}
