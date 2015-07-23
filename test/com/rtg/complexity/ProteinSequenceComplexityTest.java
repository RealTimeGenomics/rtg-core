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
package com.rtg.complexity;

import junit.framework.TestCase;

/**
 * Tests for ProteinSequenceComplexity class.
 *
 */
public class ProteinSequenceComplexityTest extends TestCase {

  public void testConstructor() {
    try {
      new ProteinSequenceComplexity(0);
      fail("Expected IAE");
    } catch (IllegalArgumentException iae) {
      // expected
    }
    try {
      new ProteinSequenceComplexity(21);
      fail("Expected IAE");
    } catch (IllegalArgumentException iae) {
      // expected
    }
    for (int i = 1; i <= 20; i++) {
      new ProteinSequenceComplexity(i);
    }
  }

  public void testComplexity() {
    final ProteinSequenceComplexity sc = new ProteinSequenceComplexity(12);
    assertEquals(0.0, sc.minComplexity("aaaaaaaaaaaaaaaaaaaaa"));
    double c = sc.minComplexity("abcdefghijklmnopqrstuvwxyz");
    assertTrue("" + c, c > 0.55 && c < 0.56);
    assertEquals(0.0, sc.maxComplexity("aaaaaaaaaaaaaaaaaaaaa"));
    c = sc.maxComplexity("abcdefghijklmnopqrstuvwxyz");
    assertTrue("" + c, c > 0.55 && c < 0.56);

    assertEquals(0.0, sc.minComplexity("aaaa"));

    c = sc.minComplexity("acgt");
    assertTrue("" + c, c > 0.26 && c < 0.27);
    c = sc.minComplexity("acgtacg");
    assertTrue("" + c, c > 0.30 && c < 0.31);
    c = sc.minComplexity("acgttttttttttttttttttttttacg");
    assertEquals(0.0, c);

    c = sc.minComplexity("acgtnxyz");
    assertTrue("" + c, c > 0.44 && c < 0.45);
  }

  public void testBadInputs() {
    final ProteinSequenceComplexity sc = new ProteinSequenceComplexity(12);
    try {
      sc.minComplexity((String) null);
      fail("Exception expected");
    } catch (NullPointerException npe) {
      // expected
    }
    try {
      sc.maxComplexity((String) null);
      fail("Exception expected");
    } catch (NullPointerException npe) {
      // expected
    }
    double c = sc.maxComplexity("");
    assertEquals(-1.0, c);
    c = sc.maxComplexity(new byte[0]);
    assertEquals(-1.0, c);
    c = sc.minComplexity("");
    assertEquals(-1.0, c);
    c = sc.minComplexity(new byte[0]);
    assertEquals(-1.0, c);
  }

  private static double complexity(String seq, int length) {
    final int[] counts = new int[27];
    for (int i = 0; i < seq.length(); i++) {
      counts[seq.charAt(i) - 'a']++;
    }

    long lf = 1;
    for (int j = 1; j <= length; j++) {
      lf *= j;
    }
    long fact = 1;
    for (int count : counts) {
      for (int k = 2; k <= count; k++) {
        fact *= k;
      }
    }
    return Math.log((double) lf / (double) fact) / Math.log(20) / length;
  }

  private static final double DELTA = 0.00001;
  public void testComplexityCalc() {
    final ProteinSequenceComplexity sc = new ProteinSequenceComplexity(12);
    final char[] seq = {'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l'};
    for (int i = 0; i < seq.length; i++) {
      seq[i] = 'a';
      final String s = new String(seq);
      final double c = sc.minComplexity(s);
      assertEquals(s, complexity(s, 12), c, DELTA);
    }
  }

  public void testCasing() {
    final ProteinSequenceComplexity sc = new ProteinSequenceComplexity(12);
    assertEquals(sc.minComplexity("acgtacgtacgtttga"), sc.minComplexity("ACGTACGTACGTTTGA"));
  }}
