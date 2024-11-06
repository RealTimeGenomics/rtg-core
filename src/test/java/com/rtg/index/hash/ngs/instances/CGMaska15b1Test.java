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
package com.rtg.index.hash.ngs.instances;

import java.io.IOException;

import com.rtg.index.hash.ngs.HashFunctionFactory;
import com.rtg.index.hash.ngs.NgsHashFunction;
import com.rtg.index.hash.ngs.ReadCall;
import com.rtg.index.hash.ngs.TemplateCall;

/**
 */
public class CGMaska15b1Test extends CGMaska1b1Test {

  @Override
  protected NgsHashFunction getHashFunction(final ReadCall readCall, final TemplateCall templateCall) {
    return new CGMaska15b1(readCall, templateCall) {
      @Override
      protected long hash(final long x) {
        //System.err.println("  hash=" + com.rtg.util.Utils.toBitsSep(x));
        return x;
      }
    };
  }

  @Override
  public void test1() throws IOException {
    final String dnar = "TGCATGCATGCATGCATGCATGCAT" +             "GCATGCATGC";
    final String dnat = "TGCATGCATGCATGCATGCATGCAT" + "aaaaaaa" + "GCATGCATGC";
    assertEquals(35, dnar.length());
    assertEquals(42, dnat.length());
    checkr(dnar, dnat, "cgmaska15test1.txt");
  }

  @Override
  public void test2() throws IOException {
    final String dnar = "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT";
    final String dnat = "TTTTTTTTTTTTTTTTTTTTTTTTT" + "aaaaaaa" + "TTTTTTTTTT";
    assertEquals(35, dnar.length());
    assertEquals(42, dnat.length());
    checkr(dnar, dnat, "cgmaska15test2.txt");
  }

  @Override
  public void test3() throws IOException {
    final String dnar = "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG";
    final String dnat = "GGGGGGGGGGGGGGGGGGGGGGGGG" + "aaaaaaa" + "GGGGGGGGGG";
    assertEquals(35, dnar.length());
    assertEquals(42, dnat.length());
    checkr(dnar, dnat, "cgmaska15test3.txt");
  }

  @Override
  public void testEmpty() throws IOException {
    checkr("", "", "cgmaska15testEmpty.txt");
  }

  @Override
  public void testFactory() {
    assertEquals(36, CGMaska15b1.FACTORY.hashBits());
    assertEquals(36, CGMaska15b1.FACTORY.windowBits());
    assertEquals(6, CGMaska15b1.FACTORY.numberWindows());
    final NgsHashFunction hf = CGMaska15b1.FACTORY.create(new ReadCallMock(null), new TemplateCallMock(null));
    assertTrue(hf != null);
    assertTrue(hf instanceof CGMaska15b1);
    assertEquals(hf.numberWindows(), CGMaska15b1.FACTORY.numberWindows());
  }

  /**
   * A factory for the outer class.
   */
  private static final HashFunctionFactory FACTORY = new CGMaska15b1.CG15HashFunctionFactory() {
    @Override
    public NgsHashFunction create(final ReadCall readCall, final TemplateCall templateCall) {
      return new CGMaska15b1(readCall, templateCall) {
        @Override
        protected long hash(long x) {
          //System.err.println("  hash=" + com.rtg.util.Utils.toBitsSep(x) + " " + x);
          return x;
        }

      };
    }
  };

  @Override
  protected HashFunctionFactory factory() {
    return FACTORY;
  }

  /** C */
  public void testCheckC() throws IOException {
    final String str =                 "acgac" + "gtgacactccgtacgta" + "ccc" +             "cgtgacaccc";
    assertEquals(35, str.length());
    //gap 7
    CheckMask.check(factory(), str,    "TcgaT" + "gtgacactccgtacgta" + "Tcc" + "AAAAAAA" + "cgtgacaccT");
    //gap 6
    CheckMask.check(factory(), str,    "TcgaT" + "gtgacactccgtacgta" + "Tcc" + "AAAAAA" +  "cgtgacaccT" + "A");
    //gap 5
    CheckMask.check(factory(), str,    "TcgaT" + "gtgacactccgtacgta" + "Tcc" + "AAAAA" +   "cgtgacaccT" + "AA");
    //gap 7 overlap 1
    final String strOv1 =              "cgacg" + "gtgacactccgtacgta" + "ccc" +             "cgtgacaccc";
    CheckMask.check(factory(), strOv1, "TcgaT" + "gtgacactccgtacgta" + "Tcc" + "AAAAAAA" + "cgtgacaccT");
    //gap 7 overlap 2
    final String strOv2 =              "gacgt" + "gtgacactccgtacgta" + "ccc" +             "cgtgacaccc";
    CheckMask.check(factory(), strOv2, "TcgaT" + "gtgacactccgtacgta" + "Tcc" + "AAAAAAA" + "cgtgacaccT");
    //gap 5 overlap 2
    CheckMask.check(factory(), strOv2, "TcgaT" + "gtgacactccgtacgta" + "Tcc" + "AAAAA" +   "cgtgacaccT" + "AA");
  }

  /** D */
  public void testCheckD() throws IOException {
    final String str =                 "acgac" + "gtg" + "acactccgtacgtaccc" +             "cgtgacaccc";
    assertEquals(35, str.length());
    //gap 7
    CheckMask.check(factory(), str,    "Tcgac" + "gtT" + "acactccgtacgtaccc" + "AAAAAAA" + "TgtgacaccT");
    //gap 6
    CheckMask.check(factory(), str,    "Tcgac" + "gtT" + "acactccgtacgtaccc" + "AAAAAA" +  "TgtgacaccT" + "A");
    //gap 5
    CheckMask.check(factory(), str,    "TcgaT" + "gtT" + "acactccgtacgtaccc" + "AAAAA" +   "TgtgacaccT" + "AA");
    //gap 7 overlap 1
    final String strOv1 =              "cgacg" + "gtg" + "acactccgtacgtaccc" +             "cgtgacaccc";
    CheckMask.check(factory(), strOv1, "Tcgac" + "gtT" + "acactccgtacgtaccc" + "AAAAAAA" + "TgtgacaccT");
    //gap 7 overlap 2
    final String strOv2 =              "gacgt" + "gtg" + "acactccgtacgtaccc" +             "cgtgacaccc";
    CheckMask.check(factory(), strOv2, "Tcgat" + "gtT" + "acactccgtacgtaccc" + "AAAAAAA" + "TgtgacaccT");
    //gap 5 overlap 2
    CheckMask.check(factory(), strOv2, "Tcgac" + "gtT" + "acactccgtacgtaccc" + "AAAAA" +   "TgtgacaccT" + "AA");
  }

  /** E */
  public void testCheckE() throws IOException {
    final String str =                 "acgac" + "gtgacact" + "ccgtacgtaccc" +             "cgtgacaccc";
    assertEquals(35, str.length());
    //gap 7
    CheckMask.check(factory(), str,    "acgac" + "TtgacacA" + "ccgtacgtaccc" + "AAAAAAA" + "TgtgacaccT");
    //gap 6
    CheckMask.check(factory(), str,    "acgac" + "TtgacacA" + "ccgtacgtaccc" + "AAAAAA" +  "TgtgacaccT" + "A");
    //gap 5
    CheckMask.check(factory(), str,    "acgac" + "TtgacacA" + "ccgtacgtaccc" + "AAAAA" +   "TgtgacaccT" + "AA");
    //gap 7 overlap 1
    final String strOv1 =              "cgacT" + "gtgacact" + "ccgtacgtaccc" +             "cgtgacaccc";
    CheckMask.check(factory(), strOv1, "acgac" + "TtgacacA" + "ccgtacgtaccc" + "AAAAAAA" + "TgtgacaccT");
    //gap 7 overlap 2
    final String strOv2 =              "gacTt" + "gtgacact" + "ccgtacgtaccc" +             "cgtgacaccc";
    CheckMask.check(factory(), strOv2, "acgac" + "TtgacacA" + "ccgtacgtaccc" + "AAAAAAA" + "TgtgacaccT");
    //gap 5 overlap 2
    CheckMask.check(factory(), strOv2, "acgac" + "TtgacacA" + "ccgtacgtaccc" + "AAAAA" +   "TgtgacaccT" + "AA");
  }

  /** F */
  public void testCheckF() throws IOException {
    final String str =                 "acgac" + "gtgacac" + "tccgtacgtaccc" +             "cgtgacaccc";
    assertEquals(35, str.length());
    //gap 7
    CheckMask.check(factory(), str,    "TcgaT" + "gtgacac" + "AccgtacgtaccT" + "AAAAAAA" + "cgtgacaccc");
    //gap 6
    CheckMask.check(factory(), str,    "acgac" + "gtgacac" + "AccgtacgtaccT" + "AAAAAA" +  "cgtgacaccc");
    //gap 5
    CheckMask.check(factory(), str,    "acgac" + "gtgacac" + "AccgtacgtaccT" + "AAAAA" +   "cgtgacaccc");
    //gap 7 overlap 1
    final String strOv1 =              "cgacg" + "gtgacac" + "tccgtacgtaccc" +             "cgtgacaccc";
    CheckMask.check(factory(), strOv1, "acgac" + "gtgacac" + "AccgtacgtaccT" + "AAAAAAA" + "cgtgacaccc");
    //gap 7 overlap 2
    final String strOv2 =              "gacgt" + "gtgacac" + "tccgtacgtaccc" +             "cgtgacaccc";
    CheckMask.check(factory(), strOv2, "acgac" + "gtgacac" + "AccgtacgtaccT" + "AAAAAAA" + "cgtgacaccc");
    //gap 5 overlap 2
    CheckMask.check(factory(), strOv2, "acgac" + "gtgacac" + "AccgtacgtaccT" + "AAAAA" +   "cgtgacaccc");
  }


}
