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
