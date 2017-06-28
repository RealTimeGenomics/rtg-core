/*
 * Copyright (c) 2017. Real Time Genomics Limited.
 *
 * Use of this source code is bound by the Real Time Genomics Limited Software Licence Agreement
 * for Academic Non-commercial Research Purposes only.
 *
 * If you did not receive a license accompanying this file, a copy must first be obtained by email
 * from support@realtimegenomics.com.  On downloading, using and/or continuing to use this source
 * code you accept the terms of that license agreement and any amendments to those terms that may
 * be made from time to time by Real Time Genomics Limited.
 */

package com.rtg.alignment;

import java.util.Arrays;
import java.util.List;

import com.rtg.util.Pair;

import junit.framework.TestCase;

/**
 * Tests the corresponding class.
 */
public class SplitAllelesTest extends TestCase {

  private String toString(final List<Pair<Integer, String[]>> split) {
    final StringBuilder sb = new StringBuilder();
    for (final Pair<Integer, String[]> s : split) {
      sb.append(s.getA()).append(' ').append(Arrays.toString(s.getB())).append('\n');
    }
    return sb.toString();
  }

  public void testNoOp() {
    final SplitAlleles sa = new SplitAlleles("");
    assertEquals("0 []\n", toString(sa.partition()));
    assertEquals(-1, sa.getColumnIndex("A"));
  }

  public void testNoAlts() {
    final SplitAlleles sa = new SplitAlleles("A");
    assertEquals("0 [A]\n", toString(sa.partition()));
    assertEquals(0, sa.getColumnIndex("A"));
  }

  public void testSimpleMatch() {
    final SplitAlleles sa = new SplitAlleles("ACGT", "ACGT");
    assertEquals("0 [ACGT, ACGT]\n", toString(sa.partition()));
    assertEquals(1, sa.getColumnIndex("ACGT"));
    assertEquals(-1, sa.getColumnIndex("A"));
  }

  public void testSimpleMatchExtra() {
    final SplitAlleles sa = new SplitAlleles("ACGT", "ACGT");
    final List<List<Pair<Integer, String[]>>> partition = sa.partition(new String[]{"TTTT"});
    assertEquals(2, partition.size());
    assertEquals("0 [ACGT, ACGT]\n", toString(partition.get(0)));
    assertEquals(1, sa.getColumnIndex("ACGT"));
    assertEquals(-1, sa.getColumnIndex("A"));
    assertEquals("0 [ACG, ACG, TTT]\n3 [T, T, T]\n", toString(partition.get(1)));
  }

  public void testSimpleMismatch() {
    final SplitAlleles sa = new SplitAlleles("AAAAA", "AACAA");
    assertEquals("0 [AA, AA]\n2 [A, C]\n3 [AA, AA]\n", toString(sa.partition()));
    assertEquals(0, sa.getColumnIndex("AAAAA"));
    assertEquals(1, sa.getColumnIndex("AACAA"));
  }

  public void testPureInsertion() {
    final SplitAlleles sa = new SplitAlleles("", "ACGT");
    assertEquals("0 [, ACGT]\n", toString(sa.partition()));
    assertEquals(0, sa.getColumnIndex(""));
  }

  public void testPureDeletion() {
    final SplitAlleles sa = new SplitAlleles("ACGT", "");
    assertEquals("0 [ACGT, ]\n", toString(sa.partition()));
  }

  public void testInsertionThenDelete() {
    // Ref  - ACGT T GC A
    // Alt0 T ACGT T GC -
    // Alt1 - ACGT C GC A
    final SplitAlleles sa = new SplitAlleles("ACGTTGCA", "TACGTTGC", "ACGTCGCA");
    assertEquals(2, sa.getColumnIndex("ACGTCGCA"));
    assertEquals("0 [, T, ]\n"
               + "0 [ACGT, ACGT, ACGT]\n"
               + "4 [T, T, C]\n"
               + "5 [GC, GC, GC]\n"
               + "7 [A, , A]\n", toString(sa.partition()));
  }

  public void testDeleteThenInsert() {
    // Ref  A CGT T GCA -
    // Alt0 - CGT T GCA T
    // Alt1 A CGT C GCA -
    final SplitAlleles sa = new SplitAlleles("ACGTTGCA", "CGTTGCAT", "ACGTCGCA");
    assertEquals("0 [A, , A]\n"
               + "1 [CGT, CGT, CGT]\n"
               + "4 [T, T, C]\n"
               + "5 [GCA, GCA, GCA]\n"
               + "8 [, T, ]\n", toString(sa.partition()));
  }

  public void testRollingDifference() {
    final SplitAlleles sa = new SplitAlleles("AAAAA", "CAAAA", "ACAAA", "AACAA", "AAACA", "AAAAC");
    assertEquals("0 [AAAAA, CAAAA, ACAAA, AACAA, AAACA, AAAAC]\n", toString(sa.partition()));
  }

  public void test1AltFromTheWild() {
    // Ref      C A GCGGAGAGACCTTTCCAGACCTGC T
    // Alt0 TAG C G GCGGAGAGACCTTTCCAGACCTGC C

    final SplitAlleles sa = new SplitAlleles("CAGCGGAGAGACCTTTCCAGACCTGCT", "TAGCGGCGGAGAGACCTTTCCAGACCTGCC");
    assertEquals("0 [, TAG]\n"
               + "0 [C, C]\n"
               + "1 [A, G]\n"
               + "2 [GCGGAGAGACCTTTCCAGACCTGC, GCGGAGAGACCTTTCCAGACCTGC]\n"
               + "26 [T, C]\n", toString(sa.partition()));
  }

  public void test2AltFromTheWild() {
    // Ref  C GCACGCG
    // Alt0 C GCACGCG GACA
    // Alt1 T GCACGCG
    final SplitAlleles sa = new SplitAlleles("CGCACGCG", "CGCACGCGGACA", "TGCACGCG");
    assertEquals("0 [C, C, T]\n"
               + "1 [GCACGCG, GCACGCG, GCACGCG]\n"
               + "8 [, GACA, ]\n", toString(sa.partition()));
  }

  public void testMatchAndDelete() {
    final SplitAlleles sa = new SplitAlleles("CGTGT", "CGT");
    assertEquals("0 [C, C]\n1 [GT, ]\n3 [GT, GT]\n", toString(sa.partition()));
  }

  public void testMatchAndInsert() {
    final SplitAlleles sa = new SplitAlleles("CGT", "CGTGT");
    assertEquals("0 [C, C]\n1 [, GT]\n1 [GT, GT]\n", toString(sa.partition()));
  }

  public void testHeterozygousInsert() {
    // Original simple trimming style gives
    // Ref  ACGT    CTGTGT
    // Alt0 ACGT TT CTGTGT
    // Alt1 ACGT CC CTGTGT

    // But doing alignment we get the alleles left aligned like this
    // Ref  ACG    T    CTGTGT
    // Alt0 ACG TT T    CTGTGT
    // Alt1 ACG    T CC CTGTGT
    final SplitAlleles sa = new SplitAlleles("ACGTCTGTGT", "ACGTTTCTGTGT", "ACGTCCCTGTGT");
    assertEquals("0 [ACG, ACG, ACG]\n"
               + "3 [, TT, ]\n"
               + "3 [T, T, T]\n"
               + "4 [, , CC]\n"
               + "4 [CTGTGT, CTGTGT, CTGTGT]\n" , toString(sa.partition()));
  }

}
