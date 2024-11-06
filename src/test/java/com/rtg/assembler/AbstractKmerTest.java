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

import com.rtg.mode.DNA;
import com.rtg.util.TestUtils;

import junit.framework.TestCase;

/**
 */
public abstract class AbstractKmerTest extends TestCase {

  protected static final String CCGT_STR = "CCGT" + "CCGT" + "CCGT" + "CCGT" + "CCGT" + "CCGT" + "CCGT" + "CCGT";
  protected static final String ACGG_STR = "ACGG" + "ACGG" + "ACGG" + "ACGG" + "ACGG" + "ACGG" + "ACGG" + "ACGG";

  protected abstract Kmer kmer(String str);

  public void test() {
    final Kmer ccgt = kmer("CCGT");
    assertEquals(4, ccgt.length());
    final Kmer acgg = kmer("ACGG");
    assertEquals("CCGT", ccgt.toString());
    assertFalse(ccgt.equals(null));
    assertFalse(ccgt.equals("false"));
    assertEquals(ccgt, acgg.reverse());
    assertTrue(0 != ccgt.hashCode());

    assertEquals('C', DNA.valueChars()[ccgt.nt(0)]);
    assertEquals('A', DNA.valueChars()[acgg.nt(0)]);
    assertEquals('G', DNA.valueChars()[ccgt.nt(2)]);
    assertEquals('G', DNA.valueChars()[acgg.nt(3)]);
    assertEquals(acgg, ccgt.minimalKmer());
    assertEquals(acgg, acgg.minimalKmer());
    assertEquals(kmer("CGGA"), acgg.successor((byte) 1));
    assertEquals(kmer("CGGT"), acgg.successor((byte) 4));
    assertEquals(kmer("CGTA"), ccgt.successor((byte) 1));
    assertEquals(kmer("CGTT"), ccgt.successor((byte) 4));

    assertEquals(kmer("AACG"), acgg.predecessor((byte) 1));
    assertEquals(kmer("TACG"), acgg.predecessor((byte) 4));
  }

  public void testEquals() {
    final Kmer empty = kmer("");
    final Kmer a = kmer("A");
    final Kmer acgt = kmer("ACGT");
    final Kmer ccgt = kmer("CCGT");
    final Kmer acgg0 = kmer("ACGG");
    final Kmer acgg1 = kmer("ACGG");
    final Kmer actt = kmer("ACTT");
    TestUtils.testOrder(new Kmer[][] {{empty}, {a}, {acgg0, acgg1}, {acgt}, {actt}, {ccgt}}, true);
  }

  public void testHashCode() {
    //regression tests
    assertEquals(-6553624, kmer("CCGT").hashCode());
    assertEquals(230748876, kmer("CCGTCCGTCCGTCCGTCCGT").hashCode());
  }

  //Test kmers of length 32 (edge case for KmerHash).
  public void test32() {
    final Kmer ccgt = kmer(CCGT_STR);
    assertEquals(32, ccgt.length());
    final Kmer acgg = kmer(ACGG_STR);
    assertEquals(CCGT_STR, ccgt.toString());
    assertFalse(ccgt.equals(null));
    assertFalse(ccgt.equals("false"));
    assertEquals(ccgt, acgg.reverse());
    assertTrue(0 != ccgt.hashCode());

    assertEquals('C', DNA.valueChars()[ccgt.nt(0)]);
    assertEquals('A', DNA.valueChars()[acgg.nt(0)]);
    assertEquals('G', DNA.valueChars()[ccgt.nt(2)]);
    assertEquals('G', DNA.valueChars()[acgg.nt(3)]);
    assertEquals('T', DNA.valueChars()[ccgt.nt(31)]);
    assertEquals('G', DNA.valueChars()[acgg.nt(31)]);
    assertEquals(acgg, ccgt.minimalKmer());
    assertEquals(acgg, acgg.minimalKmer());
    assertEquals(ACGG_STR.substring(1) + "A", acgg.successor((byte) 1).toString());
    assertEquals(ACGG_STR.substring(1) + "T", acgg.successor((byte) 4).toString());

    assertEquals("T" + ACGG_STR.substring(0, 31), acgg.predecessor((byte) 4).toString());
    assertEquals("A" + ACGG_STR.substring(0, 31), acgg.predecessor((byte) 1).toString());
  }

  public void testSuccessor() {
    checkSuccessor("ACGT");
    checkSuccessor("ACGT" + "ACGT");
    checkSuccessor(CCGT_STR);
    checkSuccessor(ACGG_STR);
  }

  protected void checkSuccessor(final String str) {
    final Kmer kmer = kmer(str);
    final String sstr = (str.length() >= 1 ? str.substring(1) : "") + "G";
    final Kmer skmer = new StringKmer(sstr);
    assertEquals(skmer, kmer.successor((byte) 3));
  }

  public void testPredecessor() {
    checkPredecessor("ACGT");
    checkPredecessor("ACGT" + "ACGT");
    checkPredecessor(CCGT_STR);
    checkPredecessor(ACGG_STR);
  }

  protected void checkPredecessor(final String str) {
    final Kmer kmer = kmer(str);
    final String sstr = "G" + (str.length() >= 1 ? str.substring(0, str.length() - 1) : "");
    final Kmer skmer = new StringKmer(sstr);
    assertEquals(str, skmer, kmer.predecessor((byte) 3));
  }

  public void testMinimal() {
    checkPredecessor("ACGT");
    checkPredecessor("ACGT" + "ACGT");
    checkPredecessor(CCGT_STR);
    checkPredecessor(ACGG_STR);
  }

  protected void checkMinimal(final String str) {
    final Kmer kmer = kmer(str);
    final Kmer krev = kmer.reverse();
    final Kmer kmin = kmer.minimalKmer();
    assertTrue(kmin.equals(krev) || kmin.equals(kmer));
  }


}
