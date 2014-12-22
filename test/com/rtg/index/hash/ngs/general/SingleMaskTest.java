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

import java.io.IOException;
import java.util.ArrayList;

import com.rtg.util.Utils;

import junit.framework.TestCase;

/**
 */
public class SingleMaskTest extends TestCase {

  public void test1() throws IOException {
    //indel 0 otherwise general case.
    final SingleMask ms = new SingleMask(4, 8, 0, 1);
    assertFalse(ms.isFrozen());
    assertEquals(8, ms.windowLength());
    assertEquals(0, ms.size());
    assertEquals(0, ms.indels());
    assertEquals(1, ms.indelLength());
    ms.integrity();
    assertEquals("Mask skeleton:liquid indels=0 indelLength=1" + LS, ms.toString());

    final Skel sk0 = new Skel(1, 2, 1);
    ms.add(sk0);
    final Skel sk1 = new Skel(7, 4, 5);
    ms.add(sk1);
    final Skel sk2 = new Skel(11, 2, 7);
    ms.add(sk2);
    ms.integrity();
    ms.freeze();
    assertTrue(ms.isFrozen());
    ms.integrity();
    assertEquals(3, ms.size());
    final String exp = ""
      + "Mask skeleton:frozen indels=0 indelLength=1" + LS
      + "0[0][1...0]<<0" + LS
      + "2[1][7...4]<<2" + LS
      + "2[2][11...10]<<6" + LS
      ;
    assertEquals(exp, ms.toString());

    assertEquals(0, ms.gaps(0));
    assertEquals(2, ms.gaps(1));
    assertEquals(2, ms.gaps(2));

    assertTrue(sk0 == ms.subSkeleton(0));
    assertTrue(sk1 == ms.subSkeleton(1));
    assertTrue(sk2 == ms.subSkeleton(2));

    final long[] res = new long[2];
    final ExtractAbstract ex = new ExtractAbstract(ms) {
      @Override
      protected void masked(final long bits) {
        res[0] = bits;
        res[1]++;
      }
    };

    ex.mask(Utils.fromBits("110011110011"), Utils.fromBits("110011110011"));
    assertEquals(1, res[1]);
    assertEquals("11111111" + "11111111", Utils.toBits(res[0]));
    res[1] = 0;

    ex.maskIndel(Utils.fromBits("110011110011"), Utils.fromBits("110011110011"));
    assertEquals(1, res[1]);
    assertEquals("11111111" + "11111111", Utils.toBits(res[0]));
    res[1] = 0;

    ex.mask(Utils.fromBits("1100001100"), Utils.fromBits("1100001100"));
    assertEquals(1, res[1]);
    assertEquals("", Utils.toBits(res[0]));
    res[1] = 0;

    ex.maskIndel(Utils.fromBits("1100001100"), Utils.fromBits("1100001100"));
    assertEquals(1, res[1]);
    assertEquals("", Utils.toBits(res[0]));
    res[1] = 0;

    ex.mask(Utils.fromBits("110110110110"), Utils.fromBits("001001001001"));
    assertEquals(1, res[1]);
    assertEquals("11101110" + "00010001", Utils.toBits(res[0]));
    res[1] = 0;

    ex.maskIndel(Utils.fromBits("110110110110"), Utils.fromBits("001001001001"));
    assertEquals(1, res[1]);
    assertEquals("11101110" + "00010001", Utils.toBits(res[0]));
    res[1] = 0;
  }

  public void test2() throws IOException {
    //indel 1 otherwise same as test1.
    final SingleMask ms = new SingleMask(4, 8, 1, 1);
    ms.integrity();
    ms.add(new Skel(1, 2, 1));
    ms.add(new Skel(7, 4, 5));
    ms.add(new Skel(11, 2, 7));
    ms.integrity();
    ms.freeze();
    ms.integrity();

    final ArrayList<Long> res = new ArrayList<>();
    final ExtractAbstract ex = new ExtractAbstract(ms) {
      @Override
      protected void masked(final long bits) {
        res.add(bits);
      }
    };

    ex.mask(Utils.fromBits("110011110011"), Utils.fromBits("110011110011"));
    assertEquals(1, res.size());
    assertEquals("11111111" + "11111111", Utils.toBits(res.get(0)));
    res.clear();

    ex.maskIndel(Utils.fromBits("110011110011"), Utils.fromBits("110011110011"));
    //unfortunately this test depends on the order the indels are generated internally
    //for (Long l : res) {
    //System.err.println(Utils.toBits(l));
    //}
    assertEquals(5, res.size());
    assertEquals("10111011" + "10111011", Utils.toBits(res.get(0)));
    assertEquals("10111111" + "10111111", Utils.toBits(res.get(1)));
    assertEquals("11111111" + "11111111", Utils.toBits(res.get(2)));
    assertEquals("1111111" + "01111111", Utils.toBits(res.get(3)));
    assertEquals("1011111" + "01011111", Utils.toBits(res.get(4)));
  }

  public void test2a() throws IOException {
    //indelLength 2 otherwise same as test1.
    final SingleMask ms = new SingleMask(4, 8, 1, 2);
    ms.integrity();
    ms.add(new Skel(1, 2, 1));
    ms.add(new Skel(7, 4, 5));
    ms.add(new Skel(11, 2, 7));
    ms.integrity();
    ms.freeze();
    ms.integrity();
    //System.err.println(ms);

    final ArrayList<Long> res = new ArrayList<>();
    final ExtractAbstract ex = new ExtractAbstract(ms) {
      @Override
      protected void masked(final long bits) {
        res.add(bits);
      }
    };

    ex.mask(Utils.fromBits("110011110011"), Utils.fromBits("110011110011"));
    assertEquals(1, res.size());
    assertEquals("11111111" + "11111111", Utils.toBits(res.get(0)));
    res.clear();

    ex.maskIndel(Utils.fromBits("110011110011"), Utils.fromBits("110011110011"));
    //unfortunately this test depends on the order the indels are generated internally
    //for (Long l : res) {
      //System.err.println(Utils.toBits(l));
    //}
    assertEquals(9, res.size());
    assertEquals("110011" + "00110011", Utils.toBits(res.get(0)));
    assertEquals("10111011" + "10111011", Utils.toBits(res.get(1)));
    assertEquals("111111" + "00111111", Utils.toBits(res.get(2)));
    assertEquals("10111111" + "10111111", Utils.toBits(res.get(3)));
    assertEquals("11111111" + "11111111", Utils.toBits(res.get(4)));
    assertEquals("1111111" + "01111111", Utils.toBits(res.get(5)));
    assertEquals("111111" + "00111111", Utils.toBits(res.get(6)));
    assertEquals("1011111" + "01011111", Utils.toBits(res.get(7)));
    assertEquals("1111" + "00001111", Utils.toBits(res.get(8)));
  }

  public void test3() throws IOException {
    //indel 3 which will hit gap boundaries - otherwise same as test1.
    final SingleMask ms = new SingleMask(4, 8, 3, 1);
    ms.integrity();
    ms.add(new Skel(1, 2, 1));
    ms.add(new Skel(7, 4, 5));
    ms.add(new Skel(11, 2, 7));
    ms.integrity();
    ms.freeze();
    ms.integrity();

    final ArrayList<Long> res = new ArrayList<>();
    final ExtractAbstract ex = new ExtractAbstract(ms) {
      @Override
      protected void masked(final long bits) {
        res.add(bits);
      }
    };

    ex.mask(Utils.fromBits("110011110011"), Utils.fromBits("110011110011"));
    assertEquals(1, res.size());
    assertEquals("11111111" + "11111111", Utils.toBits(res.get(0)));
    res.clear();

    ex.maskIndel(Utils.fromBits("110011110011"), Utils.fromBits("110011110011"));
    //unfortunately this test depends on the order the indels are generated internally
    //for (Long l : res) {
    //System.err.println(Utils.toBits(l));
    //}
    assertEquals(23, res.size());
    assertEquals("1110011" + "01110011", Utils.toBits(res.get(0)));
    assertEquals("110011" + "00110011", Utils.toBits(res.get(1)));
    assertEquals("10110011" + "10110011", Utils.toBits(res.get(2)));
    assertEquals("1111011" + "01111011", Utils.toBits(res.get(3)));
    assertEquals("111011" + "00111011", Utils.toBits(res.get(4)));

    assertEquals("10111011" + "10111011", Utils.toBits(res.get(5)));
    assertEquals("11111011" + "11111011", Utils.toBits(res.get(6)));
    assertEquals("1111011" + "01111011", Utils.toBits(res.get(7)));
    assertEquals("111111" + "00111111", Utils.toBits(res.get(8)));
    assertEquals("10111111" + "10111111", Utils.toBits(res.get(9)));

    assertEquals("11111111" + "11111111", Utils.toBits(res.get(10)));
    assertEquals("1111111" + "01111111", Utils.toBits(res.get(11)));
    assertEquals("111111" + "00111111", Utils.toBits(res.get(12)));
    assertEquals("111111" + "00111111", Utils.toBits(res.get(13)));
    assertEquals("10011111" + "10011111", Utils.toBits(res.get(14)));

    assertEquals("11011111" + "11011111", Utils.toBits(res.get(15)));
    assertEquals("1011111" + "01011111", Utils.toBits(res.get(16)));
    assertEquals("11111" + "00011111", Utils.toBits(res.get(17)));
    assertEquals("11111" + "00011111", Utils.toBits(res.get(18)));
    assertEquals("1001111" + "01001111", Utils.toBits(res.get(19)));

    assertEquals("1111" + "00001111", Utils.toBits(res.get(20)));
    assertEquals("1111" + "00001111", Utils.toBits(res.get(21)));
    assertEquals("100111" + "00100111", Utils.toBits(res.get(22)));
  }

  public void test4() throws IOException {
    //single skeleton.
    final SingleMask ms = new SingleMask(4, 2, 0, 1);
    ms.integrity();
    ms.add(new Skel(1, 2, 1));
    ms.integrity();
    ms.freeze();
    ms.integrity();
    try {
      ms.add(new Skel(1, 2, 1));
      fail();
    } catch (final RuntimeException e) {
      assertEquals("adding to frozen MaskSel", e.getMessage());
    }

    final long[] res = new long[2];
    final ExtractAbstract ex = new ExtractAbstract(ms) {
      @Override
      protected void masked(final long bits) {
        res[0] = bits;
        res[1]++;
      }
    };

    ex.mask(Utils.fromBits("10"), Utils.fromBits("01"));
    assertEquals(1, res[1]);
    assertEquals("10" + "01", Utils.toBits(res[0]));
    res[1] = 0;

    ex.maskIndel(Utils.fromBits("11"), Utils.fromBits("10"));
    assertEquals(1, res[1]);
    assertEquals("11" + "10", Utils.toBits(res[0]));
  }

  public void test5() throws IOException {
    //single skeleton - not at zero position.
    final SingleMask ms = new SingleMask(1, 2, 0, 1);
    ms.integrity();
    ms.add(new Skel(3, 2, 1));
    ms.integrity();
    ms.freeze();
    ms.integrity();

    final long[] res = new long[2];
    final ExtractAbstract ex = new ExtractAbstract(ms) {
      @Override
      protected void masked(final long bits) {
        res[0] = bits;
        res[1]++;
      }
    };

    ex.mask(Utils.fromBits("1000"), Utils.fromBits("0100"));
    assertEquals(1, res[1]);
    assertEquals("10" + "01", Utils.toBits(res[0]));
    res[1] = 0;

    ex.maskIndel(Utils.fromBits("1100"), Utils.fromBits("1000"));
    assertEquals(1, res[1]);
    assertEquals("11" + "10", Utils.toBits(res[0]));
  }
}
