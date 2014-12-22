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

import java.util.Collection;

import junit.framework.TestCase;

/**
 */
public class SkeletonTest extends TestCase {

  private void checkMasks(final Skeleton sk, final String expm) {
    final Collection<SingleMask> masks = sk.masks();
    //System.err.println(masks);
    assertEquals(sk.numberWindows(), masks.size());
    assertEquals(expm, masks.toString());
  }


  public void test1() {
    final Skeleton sk = new Skeleton(36, 12, 2, 2, 1);
    sk.integrity();
    assertEquals(36, sk.readLength());
    assertEquals(12, sk.chunkLength());
    assertEquals(2, sk.substitutions());
    assertEquals(2, sk.indels());
    assertEquals(12, sk.windowLength());
    assertEquals(24, sk.windowBits());
    final String exp = ""
      + "Inputs:" + LS
      + "  Read length   36" + LS
      + "  Window length 12" + LS
      + "  Substitutions 2" + LS
      + "  Indels        2" + LS
      + "  Indel length  1" + LS
      + "Computed:" + LS
      + "  Chunk length  12" + LS
      + "  Window chunks 1" + LS
      + "  Window actual 12" + LS
      + "  Read actual   36" + LS
      + "  Number masks  3" + LS
      ;
    assertEquals(exp, sk.toString());

    final String expm = ""
      + "[Mask skeleton:frozen indels=2 indelLength=1" + LS
      + "0[0][11...0]<<0" + LS
      + ", Mask skeleton:frozen indels=2 indelLength=1" + LS
      + "12[0][23...12]<<0" + LS
      + ", Mask skeleton:frozen indels=2 indelLength=1" + LS
      + "24[0][35...24]<<0" + LS
      + "]"
      ;
    checkMasks(sk, expm);

    final String expd = ""
      + "Mask set for readlength = 36 wordsize = 12 substitutions = 2 indels = 2 indelLength = 1" + LS
      + "" + LS
      + "000000000000000000000000111111111111" + LS
      + "000000000000111111111111000000000000" + LS
      + "111111111111000000000000000000000000" + LS
      ;
    assertEquals(expd, sk.dumpMask());
  }

  public void test1a() {
    final Skeleton sk = new Skeleton(36, 12, 2, 1, 2);
    sk.integrity();
    assertEquals(36, sk.readLength());
    assertEquals(2, sk.substitutions());
    assertEquals(1, sk.indels());
    assertEquals(2, sk.indelLength());
    assertEquals(12, sk.windowLength());
    assertEquals(24, sk.windowBits());
    final String exp = ""
      + "Inputs:" + LS
      + "  Read length   36" + LS
      + "  Window length 12" + LS
      + "  Substitutions 2" + LS
      + "  Indels        1" + LS
      + "  Indel length  2" + LS
      + "Computed:" + LS
      + "  Chunk length  12" + LS
      + "  Window chunks 1" + LS
      + "  Window actual 12" + LS
      + "  Read actual   36" + LS
      + "  Number masks  3" + LS
      ;
    assertEquals(exp, sk.toString());

    final String expm = ""
      + "[Mask skeleton:frozen indels=1 indelLength=2" + LS
      + "0[0][11...0]<<0" + LS
      + ", Mask skeleton:frozen indels=1 indelLength=2" + LS
      + "12[0][23...12]<<0" + LS
      + ", Mask skeleton:frozen indels=1 indelLength=2" + LS
      + "24[0][35...24]<<0" + LS
      + "]"
      ;
    checkMasks(sk, expm);
  }

  public void test2() {
    final Skeleton sk = new Skeleton(36, 12, 3, 2, 1);
    sk.integrity();
    assertEquals(36, sk.readLength());
    assertEquals(6, sk.chunkLength());
    assertEquals(3, sk.substitutions());
    assertEquals(2, sk.indels());
    assertEquals(12, sk.windowLength());
    assertEquals(24, sk.windowBits());
    final String exp = ""
      + "Inputs:" + LS
      + "  Read length   36" + LS
      + "  Window length 12" + LS
      + "  Substitutions 3" + LS
      + "  Indels        2" + LS
      + "  Indel length  1" + LS
      + "Computed:" + LS
      + "  Chunk length  6" + LS
      + "  Window chunks 2" + LS
      + "  Window actual 12" + LS
      + "  Read actual   30" + LS
      + "  Number masks  10" + LS
      ;
    assertEquals(exp, sk.toString());

    final String expm = ""
      + "[Mask skeleton:frozen indels=2 indelLength=1" + LS
      + "0[0][11...0]<<0" + LS
      + ", Mask skeleton:frozen indels=2 indelLength=1" + LS
      + "0[0][5...0]<<0" + LS
      + "6[1][17...12]<<6" + LS
      + ", Mask skeleton:frozen indels=2 indelLength=1" + LS
      + "0[0][5...0]<<0" + LS
      + "12[1][23...18]<<6" + LS
      + ", Mask skeleton:frozen indels=2 indelLength=1" + LS
      + "0[0][5...0]<<0" + LS
      + "18[1][29...24]<<6" + LS
      + ", Mask skeleton:frozen indels=2 indelLength=1" + LS
      + "6[0][17...6]<<0" + LS
      + ", Mask skeleton:frozen indels=2 indelLength=1" + LS
      + "6[0][11...6]<<0" + LS
      + "6[1][23...18]<<6" + LS
      + ", Mask skeleton:frozen indels=2 indelLength=1" + LS
      + "6[0][11...6]<<0" + LS
      + "12[1][29...24]<<6" + LS
      + ", Mask skeleton:frozen indels=2 indelLength=1" + LS
      + "12[0][23...12]<<0" + LS
      + ", Mask skeleton:frozen indels=2 indelLength=1" + LS
      + "12[0][17...12]<<0" + LS
      + "6[1][29...24]<<6" + LS
      + ", Mask skeleton:frozen indels=2 indelLength=1" + LS
      + "18[0][29...18]<<0" + LS
      + "]"
      ;
    checkMasks(sk, expm);
  }

  //Triggered a bug in initial testing
  public void test3() {
    final Skeleton sk = new Skeleton(3, 1, 1, 0, 1);
    sk.integrity();
    assertEquals(3, sk.readLength());
    assertEquals(1, sk.substitutions());
    assertEquals(0, sk.indels());
    assertEquals(1, sk.windowLength());
    assertEquals(2, sk.windowBits());
    final String exp = ""
      + "Inputs:" + LS
      + "  Read length   3" + LS
      + "  Window length 1" + LS
      + "  Substitutions 1" + LS
      + "  Indels        0" + LS
      + "  Indel length  1" + LS
      + "Computed:" + LS
      + "  Chunk length  1" + LS
      + "  Window chunks 1" + LS
      + "  Window actual 1" + LS
      + "  Read actual   2" + LS
      + "  Number masks  2" + LS
      ;
    assertEquals(exp, sk.toString());

    final String expm = ""
      + "[Mask skeleton:frozen indels=0 indelLength=1" + LS
      + "0[0][0...0]<<0" + LS
      + ", Mask skeleton:frozen indels=0 indelLength=1" + LS
      + "1[0][1...1]<<0" + LS
      + "]";
    checkMasks(sk, expm);
}

  //Triggered a bug in initial testing
  public void test4() {
    final Skeleton sk = new Skeleton(17, 9, 1, 0, 1);
    sk.integrity();
    assertEquals(17, sk.readLength());
    assertEquals(1, sk.substitutions());
    assertEquals(0, sk.indels());
    assertEquals(10, sk.windowLength());
    assertEquals(20, sk.windowBits());
    final String exp = ""
      + "Inputs:" + LS
      + "  Read length   17" + LS
      + "  Window length 9" + LS
      + "  Substitutions 1" + LS
      + "  Indels        0" + LS
      + "  Indel length  1" + LS
      + "Computed:" + LS
      + "  Chunk length  5" + LS
      + "  Window chunks 2" + LS
      + "  Window actual 10" + LS
      + "  Read actual   15" + LS
      + "  Number masks  3" + LS
      ;
    assertEquals(exp, sk.toString());

    final String expm = ""
      + "[Mask skeleton:frozen indels=0 indelLength=1" + LS
      + "0[0][9...0]<<0" + LS
      + ", Mask skeleton:frozen indels=0 indelLength=1" + LS
      + "0[0][4...0]<<0" + LS
      + "5[1][14...10]<<5" + LS
      + ", Mask skeleton:frozen indels=0 indelLength=1" + LS
      + "5[0][14...5]<<0" + LS
      + "]";
    checkMasks(sk, expm);
}

  //Triggered a bug in initial testing
  public void test5() {
    final Skeleton sk = new Skeleton(33, 1, 0, 0, 1);
    sk.integrity();
    assertEquals(33, sk.readLength());
    assertEquals(0, sk.substitutions());
    assertEquals(0, sk.indels());
    assertEquals(1, sk.windowLength());
    assertEquals(2, sk.windowBits());
    final String exp = ""
      + "Inputs:" + LS
      + "  Read length   33" + LS
      + "  Window length 1" + LS
      + "  Substitutions 0" + LS
      + "  Indels        0" + LS
      + "  Indel length  1" + LS
      + "Computed:" + LS
      + "  Chunk length  1" + LS
      + "  Window chunks 1" + LS
      + "  Window actual 1" + LS
      + "  Read actual   1" + LS
      + "  Number masks  1" + LS
      ;

    assertEquals(exp, sk.toString());
    final String expm = ""
      + "[Mask skeleton:frozen indels=0 indelLength=1" + LS
      + "0[0][0...0]<<0" + LS
      + "]";
    checkMasks(sk, expm);
  }

  public void testInvalid1() {
    final Skeleton sk = new Skeleton(65, 1, 0, 0, 1);
    sk.integrity();
    assertFalse(sk.valid());
    final String exp = ""
      + "Inputs:invalid" + LS
      + "  Read length   65" + LS
      + "  Window length 1" + LS
      + "  Substitutions 0" + LS
      + "  Indels        0" + LS
      + "  Indel length  1" + LS
     ;

    assertEquals(exp, sk.toString());
  }

  public void testInvalid2() {
    checkInvalid(36, 18, 3, 1, -1);
    checkInvalid(36, 18, 3, 1, 0);
    checkInvalid(-1, 1, 0, 0, 1);
    checkInvalid(65, 1, 0, 0, 1);
    checkInvalid(9, 0, 0, 0, 1);
    checkInvalid(9, 33, 0, 0, 1);
    checkInvalid(33, 33, 0, 0, 1);
    checkInvalid(35, 33, 0, 0, 1);
    checkInvalid(35, 30, 3, 0, 1);
    checkInvalid(9, 10, 0, 0, 1);
    checkInvalid(9, 4, -1, 0, 1);
    checkInvalid(9, 4, 6, 0, 1);
    checkInvalid(44, 4, 41, 0, 1);
    checkInvalid(9, 4, 2, -1, 1);
    checkInvalid(9, 4, 2, 3, 1);
  }

  private void checkInvalid(final int r, final int w, final int s, final int i, final int l) {
    final Skeleton sk = new Skeleton(r, w, s, i, l);
    sk.integrity();
    assertFalse(sk.toString(), sk.valid());
    assertTrue(sk.toString().contains(":invalid"));
  }

  public void test() {
    for (int r = 1; r <= 64; r++) {
      final int rm = Math.min(32, r);
      for (int w = 1; w < rm; w++) {
        for (int s = 0; s <= Math.min(32 - w, r - w); s++) {
          final Skeleton sk = new Skeleton(r, w, s, 0, 1);
          sk.integrity();
          //sk.masks(); //takes too long
          assertTrue(sk.valid());
        }
      }
    }

  }
}
