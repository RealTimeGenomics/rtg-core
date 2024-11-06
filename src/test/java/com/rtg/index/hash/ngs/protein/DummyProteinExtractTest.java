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
package com.rtg.index.hash.ngs.protein;


import java.io.IOException;
import java.util.ArrayList;

import com.rtg.index.hash.ngs.general.SingleMask;
import com.rtg.index.hash.ngs.general.Skel;
import com.rtg.util.Utils;

import junit.framework.TestCase;

/**
 */
public class DummyProteinExtractTest extends TestCase {

  private static class MockProteinExtract extends AbstractProteinExtract {
    protected ArrayList<Long> mResults = new ArrayList<>();
    MockProteinExtract() {
      super(new SingleMask(1, 1, 1, 1));
    }
    MockProteinExtract(SingleMask ms) {
      super(ms);
    }
    @Override
    protected void masked(long bits) {
      mResults.add(bits);
    }

  }

  public void testToString() {
    final MockProteinExtract me = new MockProteinExtract();
    assertTrue(me.integrity());
    assertEquals("ProteinExtract", me.toString());
  }

  public void test1() throws IOException {
    //indel 0 otherwise general case.
    final SingleMask ms = new SingleMask(4, 8, 0, 1);
    final Skel sk0 = new Skel(1, 2, 1);
    ms.add(sk0);
    final Skel sk1 = new Skel(7, 4, 5);
    ms.add(sk1);
    final Skel sk2 = new Skel(11, 2, 7);
    ms.add(sk2);
    ms.integrity();
    ms.freeze();
    ms.integrity();

    final MockProteinExtract ex = new MockProteinExtract(ms);

    ex.mask(new long[] {Utils.fromBits("110011110011"), Utils.fromBits("110011110011")});
    assertEquals(1, ex.mResults.size());
    assertEquals("11111111" + "11111111", Utils.toBits(ex.mResults.get(0)));
    ex.mResults.clear();

    ex.maskIndel(new long[] {Utils.fromBits("110011110011"), Utils.fromBits("110011110011"), Utils.fromBits("110011110011"), Utils.fromBits("110011110011"), Utils.fromBits("110011110011")});
    assertEquals(1, ex.mResults.size());
    assertEquals("11111111" + "11111111111111111111111111111111", Utils.toBits(ex.mResults.get(0)));
    ex.mResults.clear();

    ex.mask(new long[] {Utils.fromBits("1100001100"), Utils.fromBits("1100001100")});
    assertEquals(1, ex.mResults.size());
    assertEquals("0", Utils.toBits(ex.mResults.get(0)));
    ex.mResults.clear();

    ex.maskIndel(new long[] {Utils.fromBits("1100001100"), Utils.fromBits("1100001100"), Utils.fromBits("1100001100"), Utils.fromBits("1100001100"), Utils.fromBits("1100001100")});
    assertEquals(1, ex.mResults.size());
    assertEquals("0", Utils.toBits(ex.mResults.get(0)));
    ex.mResults.clear();

    ex.mask(new long[] {Utils.fromBits("110110110110"), Utils.fromBits("001001001001")});
    assertEquals(1, ex.mResults.size());
    assertEquals("11101110" + "00010001", Utils.toBits(ex.mResults.get(0)));
    ex.mResults.clear();

    ex.maskIndel(new long[] {Utils.fromBits("110110110110"), Utils.fromBits("001001001001"), Utils.fromBits("001001001001"), Utils.fromBits("001001001001"), Utils.fromBits("001001001001")});
    assertEquals(1, ex.mResults.size());
    assertEquals("11101110" + "00010001000100010001000100010001", Utils.toBits(ex.mResults.get(0)));
    ex.mResults.clear();
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

    final MockProteinExtract ex = new MockProteinExtract(ms);

    ex.mask(new long[] {Utils.fromBits("110011110011"), Utils.fromBits("110011110011")});
    assertEquals(1, ex.mResults.size());
    assertEquals("11111111" + "11111111", Utils.toBits(ex.mResults.get(0)));
    ex.mResults.clear();

    ex.maskIndel(new long[] {Utils.fromBits("110011110011"), Utils.fromBits("110011110011"), Utils.fromBits("110011110011"), Utils.fromBits("110011110011"), Utils.fromBits("110011110011")});
    //unfortunately this test depends on the order the indels are generated internally
    assertEquals(5, ex.mResults.size());

    assertEquals("1011101110111011101110111011101110111011", Utils.toBits(ex.mResults.get(0)));
    assertEquals("10111111" + "10111111101111111011111110111111", Utils.toBits(ex.mResults.get(1)));
    assertEquals("11111111" + "11111111111111111111111111111111", Utils.toBits(ex.mResults.get(2)));
    assertEquals("1111111" + "01111111011111110111111101111111", Utils.toBits(ex.mResults.get(3)));
    assertEquals("1011111" + "01011111010111110101111101011111", Utils.toBits(ex.mResults.get(4)));
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

    final MockProteinExtract ex = new MockProteinExtract(ms);

    ex.mask(new long[] {Utils.fromBits("110011110011"), Utils.fromBits("110011110011")});
    assertEquals(1, ex.mResults.size());
    assertEquals("11111111" + "11111111", Utils.toBits(ex.mResults.get(0)));
    ex.mResults.clear();

    ex.maskIndel(new long[] {Utils.fromBits("110011110011"), Utils.fromBits("110011110011"), Utils.fromBits("110011110011"), Utils.fromBits("110011110011"), Utils.fromBits("110011110011")});
    //unfortunately this test depends on the order the indels are generated internally
    assertEquals(9, ex.mResults.size());
    assertEquals("11001100110011001100110011001100110011", Utils.toBits(ex.mResults.get(0)));
    assertEquals("10111011" + "10111011101110111011101110111011", Utils.toBits(ex.mResults.get(1)));
    assertEquals("111111" + "00111111001111110011111100111111", Utils.toBits(ex.mResults.get(2)));
    assertEquals("10111111" + "10111111101111111011111110111111", Utils.toBits(ex.mResults.get(3)));
    assertEquals("11111111" + "11111111111111111111111111111111", Utils.toBits(ex.mResults.get(4)));
    assertEquals("1111111" + "01111111011111110111111101111111", Utils.toBits(ex.mResults.get(5)));
    assertEquals("111111" + "00111111001111110011111100111111", Utils.toBits(ex.mResults.get(6)));
    assertEquals("1011111" + "01011111010111110101111101011111", Utils.toBits(ex.mResults.get(7)));
    assertEquals("1111" + "00001111000011110000111100001111", Utils.toBits(ex.mResults.get(8)));
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

    final MockProteinExtract ex = new MockProteinExtract(ms);

    ex.mask(new long[] {Utils.fromBits("110011110011"), Utils.fromBits("110011110011")});
    assertEquals(1, ex.mResults.size());
    assertEquals("11111111" + "11111111", Utils.toBits(ex.mResults.get(0)));
    ex.mResults.clear();

    ex.maskIndel(new long[] {Utils.fromBits("110011110011"), Utils.fromBits("110011110011"), Utils.fromBits("110011110011"), Utils.fromBits("110011110011"), Utils.fromBits("110011110011")});
    //unfortunately this test depends on the order the indels are generated internally
    assertEquals(23, ex.mResults.size());
    assertEquals("111001101110011011100110111001101110011", Utils.toBits(ex.mResults.get(0)));
    assertEquals("110011" + "00110011001100110011001100110011", Utils.toBits(ex.mResults.get(1)));
    assertEquals("10110011" + "10110011101100111011001110110011", Utils.toBits(ex.mResults.get(2)));
    assertEquals("1111011" + "01111011011110110111101101111011", Utils.toBits(ex.mResults.get(3)));
    assertEquals("111011" + "00111011001110110011101100111011", Utils.toBits(ex.mResults.get(4)));

    assertEquals("10111011" + "10111011101110111011101110111011", Utils.toBits(ex.mResults.get(5)));
    assertEquals("11111011" + "11111011111110111111101111111011", Utils.toBits(ex.mResults.get(6)));
    assertEquals("1111011" + "01111011011110110111101101111011", Utils.toBits(ex.mResults.get(7)));
    assertEquals("111111" + "00111111001111110011111100111111", Utils.toBits(ex.mResults.get(8)));
    assertEquals("10111111" + "10111111101111111011111110111111", Utils.toBits(ex.mResults.get(9)));

    assertEquals("11111111" + "11111111111111111111111111111111", Utils.toBits(ex.mResults.get(10)));
    assertEquals("1111111" + "01111111011111110111111101111111", Utils.toBits(ex.mResults.get(11)));
    assertEquals("111111" + "00111111001111110011111100111111", Utils.toBits(ex.mResults.get(12)));
    assertEquals("111111" + "00111111001111110011111100111111", Utils.toBits(ex.mResults.get(13)));
    assertEquals("10011111" + "10011111100111111001111110011111", Utils.toBits(ex.mResults.get(14)));

    assertEquals("11011111" + "11011111110111111101111111011111", Utils.toBits(ex.mResults.get(15)));
    assertEquals("1011111" + "01011111010111110101111101011111", Utils.toBits(ex.mResults.get(16)));
    assertEquals("11111" + "00011111000111110001111100011111", Utils.toBits(ex.mResults.get(17)));
    assertEquals("11111" + "00011111000111110001111100011111", Utils.toBits(ex.mResults.get(18)));
    assertEquals("1001111" + "01001111010011110100111101001111", Utils.toBits(ex.mResults.get(19)));

    assertEquals("1111" + "00001111000011110000111100001111", Utils.toBits(ex.mResults.get(20)));
    assertEquals("1111" + "00001111000011110000111100001111", Utils.toBits(ex.mResults.get(21)));
    assertEquals("100111" + "00100111001001110010011100100111", Utils.toBits(ex.mResults.get(22)));
  }

  public void test4() throws IOException {
    //single skeleton.
    final SingleMask ms = new SingleMask(4, 2, 0, 1);
    ms.integrity();
    ms.add(new Skel(1, 2, 1));
    ms.integrity();
    ms.freeze();
    ms.integrity();

    final MockProteinExtract ex = new MockProteinExtract(ms);

    ex.mask(new long[] {Utils.fromBits("10"), Utils.fromBits("01")});
    assertEquals(1, ex.mResults.size());
    assertEquals("10" + "01", Utils.toBits(ex.mResults.get(0)));
    ex.mResults.clear();

    ex.maskIndel(new long[] {Utils.fromBits("11"), Utils.fromBits("10"), Utils.fromBits("00"), Utils.fromBits("00"), Utils.fromBits("00")});
    assertEquals(1, ex.mResults.size());
    assertEquals("11" + "10000000", Utils.toBits(ex.mResults.get(0)));
  }

  public void test5() throws IOException {
    //single skeleton - not at zero position.
    final SingleMask ms = new SingleMask(1, 2, 0, 1);
    ms.integrity();
    ms.add(new Skel(3, 2, 1));
    ms.integrity();
    ms.freeze();
    ms.integrity();

    final MockProteinExtract ex = new MockProteinExtract(ms);

    ex.mask(new long[] {Utils.fromBits("1000"), Utils.fromBits("0100")});
    assertEquals(1, ex.mResults.size());
    assertEquals("10" + "01", Utils.toBits(ex.mResults.get(0)));
    ex.mResults.clear();

    ex.maskIndel(new long[] {Utils.fromBits("1100"), Utils.fromBits("1000"), Utils.fromBits("0000"), Utils.fromBits("0000"), Utils.fromBits("0000")});
    assertEquals(1, ex.mResults.size());
    assertEquals("11" + "10000000", Utils.toBits(ex.mResults.get(0)));
  }
}
