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
package com.rtg.index.hash.ngs.general;

import static com.rtg.util.StringUtils.LS;

import junit.framework.TestCase;

/**
 */
public class SkeletonAltTest extends TestCase {

  public void test1() {
    final SkeletonAlt sk = new SkeletonAlt(36, 12, 2, 2);
    sk.integrity();
    final String exp = ""
      + "Inputs:" + LS
      + "  Read length   36" + LS
      + "  Window length 12" + LS
      + "  Substitutions 2" + LS
      + "  Indels        2" + LS
      + "Computed:" + LS
      + "  Chunk length  12" + LS
      + "  Small chunks  3" + LS
      + "  Big chunks    0" + LS
      + "  Window chunks 1" + LS
      + "  Window actual 12" + LS
      + "  Number masks  3" + LS
      ;

    assertEquals(exp, sk.toString());
  }

  public void test2() {
    final SkeletonAlt sk = new SkeletonAlt(36, 12, 3, 2);
    sk.integrity();
    final String exp = ""
      + "Inputs:" + LS
      + "  Read length   36" + LS
      + "  Window length 12" + LS
      + "  Substitutions 3" + LS
      + "  Indels        2" + LS
      + "Computed:" + LS
      + "  Chunk length  7" + LS
      + "  Small chunks  4" + LS
      + "  Big chunks    1" + LS
      + "  Window chunks 2" + LS
      + "  Window actual 14" + LS
      + "  Number masks  10" + LS
      ;

    assertEquals(exp, sk.toString());
  }

  //Triggered a bug in initial testing
  public void test3() {
    final SkeletonAlt sk = new SkeletonAlt(3, 1, 1, 0);
    sk.integrity();
    final String exp = ""
      + "Inputs:" + LS
      + "  Read length   3" + LS
      + "  Window length 1" + LS
      + "  Substitutions 1" + LS
      + "  Indels        0" + LS
      + "Computed:" + LS
      + "  Chunk length  1" + LS
      + "  Small chunks  1" + LS
      + "  Big chunks    1" + LS
      + "  Window chunks 1" + LS
      + "  Window actual 1" + LS
      + "  Number masks  2" + LS
      ;

    assertEquals(exp, sk.toString());
  }

  //Triggered a bug in initial testing
  public void test4() {
    final SkeletonAlt sk = new SkeletonAlt(17, 9, 1, 0);
    sk.integrity();
    final String exp = ""
      + "Inputs:" + LS
      + "  Read length   17" + LS
      + "  Window length 9" + LS
      + "  Substitutions 1" + LS
      + "  Indels        0" + LS
      + "Computed:" + LS
      + "  Chunk length  5" + LS
      + "  Small chunks  1" + LS
      + "  Big chunks    2" + LS
      + "  Window chunks 2" + LS
      + "  Window actual 11" + LS
      + "  Number masks  3" + LS
      ;

    assertEquals(exp, sk.toString());
  }

  //Triggered a bug in initial testing
  public void test5() {
    final SkeletonAlt sk = new SkeletonAlt(33, 1, 0, 0);
    sk.integrity();
    assertFalse(sk.valid());
    final String exp = ""
      + "Inputs:invalid" + LS
      + "  Read length   33" + LS
      + "  Window length 1" + LS
      + "  Substitutions 0" + LS
      + "  Indels        0" + LS
      ;

    assertEquals(exp, sk.toString());
  }

  /** Everything with a read length less than equal to 32 should be valid. */
  public void test() {
    for (int r = 1; r <= 32; r++) {
      final int rm = Math.min(32, r);
      for (int w = 1; w < rm; w++) {
        for (int s = 0; s <= Math.min(32 - w, r - w); s++) {
          final SkeletonAlt sk = new SkeletonAlt(r, w, s, 0);
          sk.integrity();
          assertTrue(sk.toString(), sk.valid());
        }
      }
    }
  }

  /** See how many cannot be satisfied above length = 32. */
  public void testInvalid() {
    int count = 0;
    int total = 0;
    for (int r = 33; r <= 64; r++) {
      final int rm = Math.min(32, r);
      for (int w = 1; w < rm; w++) {
        for (int s = 0; s <= Math.min(32 - w, r - w); s++) {
          final SkeletonAlt sk = new SkeletonAlt(r, w, s, 0);
          sk.integrity();
          total++;
          if (!sk.valid()) {
            count++;
          }
        }
      }
    }
    assertEquals(1083, count);
    assertEquals(16864, total);
  }

  public void testBest() {
    int count = 0;
    int total = 0;
    for (int r = 1; r <= 64; r++) {
      final int rm = Math.min(32, r);
      for (int w = 1; w < rm; w++) {
        for (int s = 0; s <= Math.min(32 - w, r - w); s++) {
          final SkeletonAlt ska = new SkeletonAlt(r, w, s, 0);
          final Skeleton skb = new Skeleton(r, w, s, 0, 1);
          total++;
          if (ska.valid() && ska.numberMasks() <= skb.numberWindows()) {
            count++;
          }
        }
      }
    }
    assertEquals(22816, total);
    assertEquals(21733, count);
  }

  public void testInvalid1() {
    final SkeletonAlt sk = new SkeletonAlt(65, 1, 0, 0);
    sk.integrity();
    assertFalse(sk.valid());
    final String exp = ""
      + "Inputs:invalid" + LS
      + "  Read length   65" + LS
      + "  Window length 1" + LS
      + "  Substitutions 0" + LS
      + "  Indels        0" + LS
      ;

    assertEquals(exp, sk.toString());
  }

  public void testInvalid2() {
    checkInvalid(-1, 1, 0, 0);
    checkInvalid(65, 1, 0, 0);
    checkInvalid(65, 31, 0, 0);
    checkInvalid(9, 0, 0, 0);
    checkInvalid(9, 33, 0, 0);
    checkInvalid(33, 33, 0, 0);
    checkInvalid(35, 33, 0, 0);
    checkInvalid(35, 30, 3, 0);
    checkInvalid(9, 10, 0, 0);
    checkInvalid(9, 4, -1, 0);
    checkInvalid(9, 4, 6, 0);
    checkInvalid(44, 4, 41, 0);
    checkInvalid(9, 4, 2, -1);
    checkInvalid(9, 4, 2, 3);
  }

  public void testInvalid3() {
    checkInvalid(65, 1, 0, 0);
    checkInvalid(65, 31, 0, 0);
  }

  private void checkInvalid(final int r, final int w, final int s, final int i) {
    final SkeletonAlt sk = new SkeletonAlt(r, w, s, i);
    sk.integrity();
    assertFalse(sk.toString(), sk.valid());
    assertTrue(sk.toString().contains(":invalid"));
  }

}
