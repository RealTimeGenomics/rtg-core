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
package com.rtg.index.hash.ngs.instances;

import java.io.IOException;

import com.rtg.index.hash.ngs.HashFunctionFactory;
import com.rtg.index.hash.ngs.NgsHashFunction;
import com.rtg.index.hash.ngs.ReadCall;
import com.rtg.index.hash.ngs.TemplateCall;

/**
 */
public class CGMaska15b1altTest extends CGMaska1b1altTest {

  @Override
  protected NgsHashFunction getHashFunction(final ReadCall readCall, final TemplateCall templateCall) {
    return new CGMaska15b1alt(readCall, templateCall) {
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
    //System.err.println(displayDNA(dnat));
    assertEquals(35, dnar.length());
    assertEquals(42, dnat.length());
    checkr(dnar, dnat, "cgmaska15alttest1.txt");
  }

  @Override
  public void test2() throws IOException {
    final String dnar = "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT";
    final String dnat = "TTTTTTTTTTTTTTTTTTTTTTTTT" + "aaaaaaa" + "TTTTTTTTTT";
    //System.err.println(displayDNA(dnat));
    assertEquals(35, dnar.length());
    assertEquals(42, dnat.length());
    checkr(dnar, dnat, "cgmaska15test2.txt");
  }

  @Override
  public void test3() throws IOException {
    final String dnar = "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG";
    final String dnat = "GGGGGGGGGGGGGGGGGGGGGGGGG" + "aaaaaaa" + "GGGGGGGGGG";
    //System.err.println(displayDNA(dnat));
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
    assertEquals(36, CGMaska15b1alt.FACTORY.hashBits());
    assertEquals(36, CGMaska15b1alt.FACTORY.windowBits());
    assertEquals(6, CGMaska15b1alt.FACTORY.numberWindows());
    final NgsHashFunction hf = CGMaska15b1alt.FACTORY.create(new ReadCallMock(null), new TemplateCallMock(null));
    assertTrue(hf != null);
    assertTrue(hf instanceof CGMaska15b1alt);
    assertEquals(hf.numberWindows(), CGMaska15b1alt.FACTORY.numberWindows());
  }

  /**
   * A factory for the outer class.
   */
  private static final HashFunctionFactory FACTORY = new CGMaska15b1alt.CG15HashFunctionFactory() {
    @Override
    public NgsHashFunction create(final ReadCall readCall, final TemplateCall templateCall) {
      return new CGMaska15b1alt(readCall, templateCall) {
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
    //added Alt suffix for suppress cpd
    final String strAlt =                 "acgac" + "gtgacactccgtacgta" + "ccc" +             "cgtgacaccc";
    assertEquals(35, strAlt.length());
    //gap 7
    CheckMask.check(factory(), strAlt,    "TcgaT" + "gtgacactccgtacgta" + "Tcc" + "AAAAAAA" + "cgtgacaccT");
    //gap 6
    CheckMask.check(factory(), strAlt,    "TcgaT" + "gtgacactccgtacgta" + "Tcc" + "AAAAAA" +  "cgtgacaccT" + "A");
    //gap 5
    CheckMask.check(factory(), strAlt,    "TcgaT" + "gtgacactccgtacgta" + "Tcc" + "AAAAA" +   "cgtgacaccT" + "AA");
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
    final String strAlt =                 "acgac" + "gtg" + "acactccgtacgtaccc" +             "cgtgacaccc";
    assertEquals(35, strAlt.length());
    //gap 7
    CheckMask.check(factory(), strAlt,    "Tcgac" + "gtT" + "acactccgtacgtaccc" + "AAAAAAA" + "TgtgacaccT");
    //gap 6
    CheckMask.check(factory(), strAlt,    "Tcgac" + "gtT" + "acactccgtacgtaccc" + "AAAAAA" +  "TgtgacaccT" + "A");
    //gap 5
    CheckMask.check(factory(), strAlt,    "TcgaT" + "gtT" + "acactccgtacgtaccc" + "AAAAA" +   "TgtgacaccT" + "AA");
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
    //CheckMask.check(factory(), str,    "acgac" + "TtgacacA" + "ccgtacgtaccc" + "AAAAAAA" + "TgtgacaccT");
    //gap 6
    //CheckMask.check(factory(), str,    "acgac" + "TtgacacA" + "ccgtacgtaccc" + "AAAAAA" +  "TgtgacaccT" + "A");
    //gap 5
    //CheckMask.check(factory(), str,    "acgac" + "TtgacacA" + "ccgtacgtaccc" + "AAAAA" +   "TgtgacaccT" + "AA");

    //gap 7 overlap 1
    final String strOv1 =              "cgacT" + "gtgacact" + "ccgtacgtaccc" +             "cgtgacaccc";
    CheckMask.check(factory(), strOv1, "acgac" + "TtgacacA" + "ccgtacgtaccc" + "AAAAAAA" + "TgtgacaccT");

    final String strOv2 =              "gacTt" + "gtgacact" + "ccgtacgtaccc" +             "cgtgacaccc";
    //gap 7 overlap 2
    CheckMask.check(factory(), strOv2, "acgac" + "TtgacacA" + "ccgtacgtaccc" + "AAAAAAA" + "TgtgacaccT");
    //gap 5 overlap 2
    CheckMask.check(factory(), strOv2, "acgac" + "TtgacacA" + "ccgtacgtaccc" + "AAAAA" +   "TgtgacaccT" + "AA");

    final String strOv3 =              "acTtg" + "gtgacact" + "ccgtacgtaccc" +             "cgtgacaccc";
    //gap 7 overlap 3
    CheckMask.check(factory(), strOv3,    "acgac" + "TtgacacA" + "ccgtacgtaccc" + "AAAAAAA" + "TgtgacaccT");
    //gap 6 overlap 3
    CheckMask.check(factory(), strOv3,    "acgac" + "TtgacacA" + "ccgtacgtaccc" + "AAAAAA" +  "TgtgacaccT" + "A");
    //gap 5 overlap 3
    CheckMask.check(factory(), strOv3,    "acgac" + "TtgacacA" + "ccgtacgtaccc" + "AAAAA" +   "TgtgacaccT" + "AA");
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
