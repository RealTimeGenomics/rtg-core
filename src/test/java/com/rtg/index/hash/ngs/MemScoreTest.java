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
package com.rtg.index.hash.ngs;

import java.util.Arrays;

import com.rtg.mode.DNA;
import com.rtg.util.PortableRandom;

import junit.framework.TestCase;

/**
 * Tests the corresponding class.
 */
public class MemScoreTest extends TestCase {

  public void test() {
    for (int k = 1; k < 200; k += 17) {
      final byte[] seq = new byte[k];
      assertEquals(0, MemScore.score(seq, 1, seq, 1, false));
      assertEquals(0, MemScore.score(seq, 1, seq, 1, false));
      assertEquals(0, MemScore.score(seq, seq.length, seq, 1, true));
      if (k > 5) {
        assertEquals(0, MemScore.score(seq, 1, seq, 1, false));
        assertEquals(5, MemScore.score(seq, 6, seq, 1, false));
        final byte[] seq1 = new byte[k];
        seq1[3] = 1;
        assertEquals(1, MemScore.score(seq1, 1, seq, 1, false));
      }
      Arrays.fill(seq, (byte) 1);
      assertEquals(0, MemScore.score(seq, 1, seq, 1, false));
      assertEquals(k, MemScore.score(seq, seq.length, seq, 1, true));
    }
  }

  private static byte[] getRandom(final PortableRandom r, final int length) {
    final byte[] template = new byte[length];
    r.nextBytes(template);
    for (int k = 0; k < template.length; ++k) {
      template[k] &= 3;
    }
    return template;
  }

  public void testRandom() {
    final PortableRandom r = new PortableRandom();
    final byte[] template = getRandom(r, 10000);
    for (int k = 1; k < 10000; k += k) {
      final byte[] read = getRandom(r, k);
      final int tDelta = r.nextInt(template.length);
      final int rDelta = r.nextInt(read.length);
      final int score = MemScore.score(template, 1 + tDelta, read, 1 + rDelta, false);
      assertTrue(score >= 0);
      assertTrue(score <= read.length);
      final byte[] rc = new byte[read.length];
      for (int j = 0, i = read.length - 1; j < read.length; ++j, --i) {
        rc[j] = DNA.complement(read[i]);
      }
      assertEquals(score, MemScore.score(template, 1 + tDelta, rc, read.length - rDelta , true));
    }
  }

  public void testCycle() {
    final byte[] seq = new byte[50];
    for (int k = 0; k < seq.length; ++k) {
      seq[k] = (byte) (1 + (k & 3));
    }
    assertEquals(0, MemScore.score(seq, 1, seq, 1, false));
    final byte[] seq1 = new byte[25];
    for (int k = 0; k < seq1.length; ++k) {
      seq1[k] = (byte) (1 + (k & 3));
    }
    assertEquals(0, MemScore.score(seq, 1, seq1, 1, false));
    assertEquals(25, MemScore.score(seq, 2, seq1, 1, false));
    assertEquals(25, MemScore.score(seq, 3, seq1, 1, false));
    assertEquals(25, MemScore.score(seq, 4, seq1, 1, false));
    assertEquals(0, MemScore.score(seq, 5, seq1, 1, false));
  }

  private void check(final long a0, final long a1, final long b0, final long b1, final int expected) {
    final long diff = MemScore.simpleDiff(a0, a1, b0, b1);
    final int res = Long.bitCount(diff);
    assertEquals("a0=" + a0 + " a1=" + a1 + " b0=" + b0 + " b1=" + b1, expected, res);
  }

  public void testSimpleScore() {
    check(0, 0, 0, 0, 0);
    check(0, 0, 0, 1, 1);
    check(0, 0, 1, 0, 1);
    check(0, 0, 1, 1, 1);
    check(0, 1, 0, 0, 1);
    check(0, 1, 0, 1, 0);
    check(0, 1, 1, 0, 1);
    check(0, 1, 1, 1, 1);
    check(1, 0, 0, 0, 1);
    check(1, 0, 0, 1, 1);
    check(1, 0, 1, 0, 0);
    check(1, 0, 1, 1, 1);
    check(1, 1, 0, 0, 1);
    check(1, 1, 0, 1, 1);
    check(1, 1, 1, 0, 1);
    check(1, 1, 1, 1, 0);
    check(-1L, 0, 0, 0, 64);
    check(0, -1L, 0, 0, 64);
    check(0, 0, -1L, 0, 64);
    check(0, 0, 0, -1L, 64);
    check(0, -1L, 0, -1L, 0);
  }

  private void checkExact(final MemScore ms, final long a0, final long a1, final long b0, final long b1, final int expected) {
    final int res = ms.fastScore(a0, a1, b0, b1);
    assertEquals("a0=" + a0 + " a1=" + a1 + " b0=" + b0 + " b1=" + b1, expected, res);
    final int indel = ms.indelScore(a0, a1, b0, b1);
    assertTrue("a0=" + a0 + " a1=" + a1 + " b0=" + b0 + " b1=" + b1 + " res=" + res + " indel=" + indel, indel <= res);
  }

  public void check(final int l, final long mask) {
    final MemScore ms = new MemScore(l);
    ms.integrity();
    checkExact(ms, 0, 0, 0, 0, 0);
    checkExact(ms, 0, 0, 0, 1, 1);
    checkExact(ms, 0, 0, 1, 0, 1);
    checkExact(ms, 0, 0, 1, 1, 1);
    checkExact(ms, 0, 1, 0, 0, 1);
    checkExact(ms, 0, 1, 0, 1, 0);
    checkExact(ms, 0, 1, 1, 0, 1);
    checkExact(ms, 0, 1, 1, 1, 1);
    checkExact(ms, 1, 0, 0, 0, 1);
    checkExact(ms, 1, 0, 0, 1, 1);
    checkExact(ms, 1, 0, 1, 0, 0);
    checkExact(ms, 1, 0, 1, 1, 1);
    checkExact(ms, 1, 1, 0, 0, 1);
    checkExact(ms, 1, 1, 0, 1, 1);
    checkExact(ms, 1, 1, 1, 0, 1);
    checkExact(ms, 1, 1, 1, 1, 0);
    checkExact(ms, mask, 0, 0, 0, l);
    checkExact(ms, 0, mask, 0, 0, l);
    checkExact(ms, 0, 0, -1L, 0, l);
    checkExact(ms, 0, 0, 0, -1L, l);
    checkExact(ms, 0, mask, 0, -1L, 0);
  }

  public void testAll() {
    for (int i = 1; i < 64; ++i) {
      check(i, (1L << i) - 1);
    }
    check(64, -1L);
  }

  public void test2() {
    final MemScore ms = new MemScore(36);
    assertEquals(2, ms.fastScore(0L, 0L, 2, 1));
    assertEquals(1, ms.fastScore(2L, 0L, 2, 1));
    assertEquals("Score length=36", ms.toString());
  }

  private void checkIndel(final MemScore ms, final long a0, final long a1, final long b0, final long b1, final int exact, final int indelex) {
    final int res = ms.fastScore(a0, a1, b0, b1);
    final int indel = ms.indelScore(a0, a1, b0, b1);
    assertEquals("a0=" + a0 + " a1=" + a1 + " b0=" + b0 + " b1=" + b1, exact, res);
    assertEquals("a0=" + a0 + " a1=" + a1 + " b0=" + b0 + " b1=" + b1, indelex, indel);
  }

  private void checkIndel(final MemScore ms, final String a0, final String a1, final String b0, final String b1, final int exact, final int indelex) {
    checkIndel(ms, Long.valueOf(a0, 2), Long.valueOf(a1, 2), Long.valueOf(b0, 2), Long.valueOf(b1, 2), exact, indelex);
    checkIndel(ms, Long.valueOf(b0, 2), Long.valueOf(b1, 2), Long.valueOf(a0, 2), Long.valueOf(a1, 2), exact, indelex);
  }

  /**
   * Template t[actgac]
   *          1 001100
   *          1 011001
   */
  public void testIndel6() {
    final MemScore ms = new MemScore(6);
    //Exact actgac
    //      001100
    //      011001
    checkIndel(ms, "001100", "011001", "001100", "011001", 0, 0);

    //Deletion tactac
    //         100100
    //         101101
    checkIndel(ms, "001100", "011001", "100100", "101101", 4, 2);

    //Insertion tactggac
    //          10011100
    //          10110001
    checkIndel(ms, "001100", "011001", "10011100", "10110001", 3, 2);
  }

  /**
   * Template t[actgac]
   *          1 001100
   *          1 011001
   */
  public void testIndel6s() {
    final MemScore ms = new MemScore(6);
    //Subs  actgag
    //      001101
    //      011000
    checkIndel(ms, "001100", "011001", "001101", "011000", 1, 1);

    //Deletion tactag
    //         100101
    //         101100
    checkIndel(ms, "001100", "011001", "100101", "101100", 5, 3);

    //Insertion tactggac
    //          10011101
    //          10110000
    checkIndel(ms, "001100", "011001", "10011101", "10110000", 4, 2);
  }

  /**
   * Template a[actgac]
   *          0 001100
   *          0 011001
   */
  public void testIndel() {
    for (int i = 6; i < 64; ++i) {
      final MemScore ms = new MemScore(i);
      //Exact actgac
      //      001100
      //      011001
      checkIndel(ms, "001100", "011001", "001100", "011001", 0, 0);

      //Deletion aactac
      //         000100
      //         001101
      checkIndel(ms, "001100", "011001", "000100", "001101", 3, 2);

      //Insertion aactggac
      //          00011100
      //          00110001
      checkIndel(ms, "001100", "011001", "00011100", "00110001", 3, 2);
    }
  }

  public void test63() {
    final MemScore ms = new MemScore(63);
    checkExact(ms, -1L >>> 1, -1L >>> 1, -1L, -1L, 0);
  }

  public void testIndel63() {
    final MemScore ms = new MemScore(63);
    checkIndel(ms, -1L >>> 1, -1L >>> 1, -1L, -1L, 0, 0);
    checkIndel(ms, -1L >>> 2, -1L >>> 1, -1L, -1L, 1, 1);
  }

  /**
   * Bug found in the wild.
   *     ACCGTGCGTAATTTTTTATCACGGCTTaTACTTCAT
   * cgttACCGTGCGTAATTTTTTATCACGGCTTTACTTCATcc
   */
  public void testBug36() {
    final MemScore ms = new MemScore(36);
    checkIndel(
        ms,
            "000111011001111110100011011010011001", //r0
            "011010101001111110110100111010111101", //r1

        "0111000111011001111110100011011100110010", //t0
        "1011011010101001111110110100111101111011", //t1
        8,
        2
    );
  }
}

