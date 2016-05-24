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

import static com.rtg.util.StringUtils.LS;

import java.io.IOException;

import com.rtg.index.hash.ngs.HashFunctionFactory;
import com.rtg.index.hash.ngs.NgsHashFunction;
import com.rtg.index.hash.ngs.ReadCall;
import com.rtg.index.hash.ngs.TemplateCall;

import org.junit.Assert;

/**
 */
public class CGMaska1b1Test extends AbstractSplitTest {

  @Override
  protected NgsHashFunction getHashFunction(final ReadCall readCall, final TemplateCall templateCall) {
    return new CGMaska1b1(readCall, templateCall) {
      @Override
      protected long hash(final long x) {
        //System.err.println("  hash=" + com.rtg.util.Utils.toBitsSep(x));
        return x;
      }
    };
  }

  private static final String PREFIX = ""
    + "CG l=35 w=18 e=1" + LS
    + "number windows=2" + LS
    + "set name=1234 length=23" + LS;

  private static final String SUFFIX = "done" + LS;

  private static final String EXPECTED_1 = ""
    + PREFIX
    + "readCall index=0 id=5" + LS
    + " hash=00000000:00000000:00000000:00001100:11001100:11001110:10101010:10101010" + LS
    + "readCall index=1 id=5" + LS
    + " hash=00000000:00000000:00000000:00000001:10011001:10011001:01010101:01010101" + LS
    + "templateCall position=2 index=0" + LS
    + " hash=00000000:00000000:00000000:00001100:11001100:11001110:10101010:10101010" + LS
    + "templateCall position=2 index=0" + LS
    + " hash=00000000:00000000:00000000:00001001:11001100:11001101:01001010:10101010" + LS
    + "templateCall position=2 index=0" + LS
    + " hash=00000000:00000000:00000000:00000011:01001100:11001110:10101010:10101010" + LS
    + "templateCall position=2 index=1" + LS
    + " hash=00000000:00000000:00000000:00000001:10011001:10011001:01010101:01010101" + LS
    + "templateCall position=2 index=1" + LS
    + " hash=00000000:00000000:00000000:00000011:00101001:10011000:10101001:01010101" + LS
    + "templateCall position=2 index=1" + LS
    + " hash=00000000:00000000:00000000:00000110:01001001:10011001:01010001:01010101" + LS
    + SUFFIX
    ;
  public void test1() throws IOException {
    final String dnar = "TGCATGCATGCATGCATGCATGCAT" +             "GCATGCATGC";
    final String dnat = "TGCATGCATGCATGCATGCATGCAT" + "aaaaaaa" + "GCATGCATGC";
    assertEquals(35, dnar.length());
    assertEquals(42, dnat.length());
    check(dnar, dnat, EXPECTED_1);
  }

  private static final String EXPECTED_2 = ""
    + PREFIX
    + "readCall index=0 id=5" + LS
    + " hash=00000000:00000000:00000000:00001111:11111111:11111111:11111111:11111111" + LS
    + "readCall index=1 id=5" + LS
    + " hash=00000000:00000000:00000000:00000111:11111111:11111101:11111111:11111111" + LS
    + "templateCall position=2 index=0" + LS
    + " hash=00000000:00000000:00000000:00001111:11111111:11111111:11111111:11111111" + LS
    + "templateCall position=2 index=0" + LS
    + " hash=00000000:00000000:00000000:00001111:11111111:11111111:11111111:11111111" + LS
    + "templateCall position=2 index=0" + LS
    + " hash=00000000:00000000:00000000:00001111:11111111:11111111:11111111:11111111" + LS
    + "templateCall position=2 index=1" + LS
    + " hash=00000000:00000000:00000000:00000111:11111111:11111101:11111111:11111111" + LS
    + "templateCall position=2 index=1" + LS
    + " hash=00000000:00000000:00000000:00000111:11101111:11111101:11111011:11111111" + LS
    + "templateCall position=2 index=1" + LS
    + " hash=00000000:00000000:00000000:00000111:11001111:11111101:11110011:11111111" + LS
    + SUFFIX
    ;
  public void test2() throws IOException {
    final String dnar = "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT";
    final String dnat = "TTTTTTTTTTTTTTTTTTTTTTTTT" + "aaaaaaa" + "TTTTTTTTTT";
    assertEquals(35, dnar.length());
    assertEquals(42, dnat.length());
    check(dnar, dnat, EXPECTED_2);
  }

  private static final String EXPECTED_3 = ""
    + PREFIX
    + "readCall index=0 id=5" + LS
    + " hash=00000000:00000000:00000000:00001111:11111111:11111100:00000000:00000000" + LS
    + "readCall index=1 id=5" + LS
    + " hash=00000000:00000000:00000000:00000111:11111111:11111100:00000000:00000000" + LS

    + "templateCall position=2 index=0" + LS
    + " hash=00000000:00000000:00000000:00001111:11111111:11111100:00000000:00000000" + LS
    + "templateCall position=2 index=0" + LS
    + " hash=00000000:00000000:00000000:00001111:11111111:11111100:00000000:00000000" + LS
    + "templateCall position=2 index=0" + LS
    + " hash=00000000:00000000:00000000:00001111:11111111:11111100:00000000:00000000" + LS

    + "templateCall position=2 index=1" + LS
    + " hash=00000000:00000000:00000000:00000111:11111111:11111100:00000000:00000000" + LS
    + "templateCall position=2 index=1" + LS
    + " hash=00000000:00000000:00000000:00000111:11101111:11111100:00000000:00000000" + LS
    + "templateCall position=2 index=1" + LS
    + " hash=00000000:00000000:00000000:00000111:11001111:11111100:00000000:00000000" + LS
    + SUFFIX
    ;
  public void test3() throws IOException {
    final String dnar = "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG";
    final String dnat = "GGGGGGGGGGGGGGGGGGGGGGGGG" + "aaaaaaa" + "GGGGGGGGGG";
    assertEquals(35, dnar.length());
    assertEquals(42, dnat.length());
    check(dnar, dnat, EXPECTED_3);
  }

  public void testEmpty() throws IOException {
    check("", "", PREFIX);
  }

  public void testFactory() {
    assertEquals(36, CGMaska1b1.FACTORY.hashBits());
    assertEquals(36, CGMaska1b1.FACTORY.windowBits());
    assertEquals(2, CGMaska1b1.FACTORY.numberWindows());
    final NgsHashFunction hf = CGMaska1b1.FACTORY.create(new ReadCallMock(null), new TemplateCallMock(null));
    assertTrue(hf != null);
    assertTrue(hf instanceof CGMaska1b1);
    assertEquals(hf.numberWindows(), CGMaska1b1.FACTORY.numberWindows());
  }

  /**
   * A factory for the outer class.
   */
  private static final HashFunctionFactory FACTORY = new CGMaska1b1.CG1HashFunctionFactory() {
    @Override
    public NgsHashFunction create(final ReadCall readCall, final TemplateCall templateCall) {
      return new CGMaska1b1(readCall, templateCall) {
        @Override
        protected long hash(long x) {
          //System.err.println("  hash=" + com.rtg.util.Utils.toBitsSep(x) + " " + x);
          return x;
        }

      };
    }
  };

  protected HashFunctionFactory factory() {
    return FACTORY;
  }

  /**
   * Check that all 0, 1 substitutions on the string are found.
   */
  public void testAllSubstitutions7() throws IOException {
    final String str = "acgacctgacacccgtacgtaccccgtgacacccg";
    assertEquals(35, str.length());
    final SubstituteCG sub = new SubstituteCG(str, factory(), true, 7, 1);
    final int missed = sub.substituteProtected(1);
    Assert.assertEquals(sub.toString(), 0, missed);
  }

  /**
   * Check that all 0, 1, 2 substitutions on the string are found.
   */
  public void testAllSubstitutions6() throws IOException {
    final String str = "acgacctgacacccgtacgtaccccgtgacacccg";
    assertEquals(35, str.length());
    final SubstituteCG sub = new SubstituteCG(str, factory(), true, 6, 1);
    Assert.assertEquals(0, sub.substituteProtected(1));
  }

  /**
   * Check that all 0, 1, 2 substitutions on the string are found.
   */
  public void testAllSubstitutions5() throws IOException {
    final String str = "acgacctgacacccgtacgtaccccgtgacacccg";
    assertEquals(35, str.length());
    final SubstituteCG sub = new SubstituteCG(str, factory(), true, 5, 1);
    Assert.assertEquals(0, sub.substituteProtected(1));
  }

  /**
   * Check that all 0, 1, 2 substitutions on the string are found.
   * @throws IOException
   */
  public void testIndel() throws IOException {
    final String str = "acgacgtgacacccgtacgtaccccgtgacacccg";
    assertEquals(35, str.length());
    SubstituteIndel.checkIndel(factory(), str, 1, 7);
  }

  /** Exact */
  public void testCheckExact() throws IOException {
    final String str = "acgac" + "gtgacactccgtacgtaccc" + "cgtgacaccc";
    assertEquals(35, str.length());
    //gap 7 overlap 0
    CheckMask.check(factory(), str, "acgac" + "gtgacactccgtacgtaccc" + "AAAAAAA" + "cgtgacaccc");
    //gap 6 overlap 0
    CheckMask.check(factory(), str, "acgac" + "gtgacactccgtacgtaccc" + "AAAAAA" + "cgtgacaccc" + "A");
    //gap 5 overlap 0
    CheckMask.check(factory(), str, "acgac" + "gtgacactccgtacgtaccc" + "AAAAA" + "cgtgacaccc" + "AA");
    //gap 7 overlap 1
    final String strOv1 = "cgacg" + "gtgacactccgta" + "cgtaccc" + "cgtgacaccc";
    CheckMask.check(factory(), strOv1, "acgac" + "gtgacactccgtacgtaccc" + "AAAAAAA" + "cgtgacaccc");
    //gap 7 overlap 2
    final String strOv2 = "gacgt" + "gtgacactccgta" + "cgtaccc" + "cgtgacaccc";
    CheckMask.check(factory(), strOv2, "acgac" + "gtgacactccgtacgtaccc" + "AAAAAAA" + "cgtgacaccc");
    //gap 5 overlap 2
    CheckMask.check(factory(), strOv2, "acgac" + "gtgacactccgtacgtaccc" + "AAAAA" + "cgtgacaccc" + "AA");
  }

  /** A */
  public void testCheckA() throws IOException {
    final String str =                 "acgac" + "gtgacactccgta" + "cgtaccc" +             "cgtgacaccc";
    assertEquals(35, str.length());
    //gap 7 overlap 0
    CheckMask.check(factory(), str,    "acgac" + "gtgacactccgta" + "Tgtaccc" + "AAAAAAA" + "cgtgacaccT");
    //gap 6 overlap 0
    CheckMask.check(factory(), str,    "acgac" + "gtgacactccgta" + "Tgtaccc" + "AAAAAA" +  "cgtgacaccT" + "A");
    //gap 5 overlap 0
    CheckMask.check(factory(), str,    "acgac" + "gtgacactccgta" + "Tgtaccc" + "AAAAA" +   "cgtgacaccT" + "AA");
    //gap 7 overlap 1
    final String strOv1 =              "cgacg" + "gtgacactccgta" + "cgtaccc" +             "cgtgacaccc";
    CheckMask.check(factory(), strOv1, "acgac" + "gtgacactccgta" + "Tgtaccc" + "AAAAAAA" + "cgtgacaccT");
    //gap 7 overlap 2
    final String strOv2 =              "gacgt" + "gtgacactccgta" + "cgtaccc" +             "cgtgacaccc";
    CheckMask.check(factory(), strOv2, "acgac" + "gtgacactccgta" + "Tgtaccc" + "AAAAAAA" + "cgtgacaccT");
    //gap 5 overlap 2
    CheckMask.check(factory(), strOv2, "acgac" + "gtgacactccgta" + "Tgtaccc" + "AAAAA" +   "cgtgacaccT" + "AA");
  }

  /** B */
  public final void testCheckB() throws IOException {
    final String str =                 "acgac" + "gtgacactccgta" + "cgtaccc" +             "cgtgacaccc";
    assertEquals(35, str.length());
    //gap 7
    CheckMask.check(factory(), str,    "Tcgac" + "gtgacactccgtT" + "cgtaccc" + "AAAAAAA" + "cgtgacaccc");
    //gap 6
    CheckMask.check(factory(), str,    "Tcgac" + "gtgacactccgtT" + "cgtaccc" + "AAAAAA" +  "cgtgacaccc");
    //gap 5
    CheckMask.check(factory(), str,    "Tcgac" + "gtgacactccgtT" + "cgtaccc" + "AAAAA" +   "cgtgacaccc");
    //gap 7 overlap 1
    final String strOv1 =              "cgacg" + "gtgacactccgta" + "cgtaccc" +             "cgtgacaccc";
    CheckMask.check(factory(), strOv1, "Tcgac" + "gtgacactccgtT" + "cgtaccc" + "AAAAAAA" + "cgtgacaccc");
    //gap 7 overlap 2
    final String strOv2 =              "gacgt" + "gtgacactccgta" + "cgtaccc" +             "cgtgacaccc";
    CheckMask.check(factory(), strOv2, "Tcgac" + "gtgacactccgtT" + "cgtaccc" + "AAAAAAA" + "cgtgacaccc");
    //gap 5 overlap 2
    CheckMask.check(factory(), strOv2, "Tcgac" + "gtgacactccgtT" + "cgtaccc" + "AAAAA" +   "cgtgacaccc");
  }

}
