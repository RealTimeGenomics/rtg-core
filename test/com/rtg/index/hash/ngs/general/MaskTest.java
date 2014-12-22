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
import java.io.StringWriter;

import com.rtg.index.hash.ngs.HashFunctionFactory;
import com.rtg.index.hash.ngs.NgsHashFunction;
import com.rtg.index.hash.ngs.ReadCall;
import com.rtg.index.hash.ngs.TemplateCall;
import com.rtg.index.hash.ngs.instances.AbstractSplitTest;
import com.rtg.index.hash.ngs.instances.Substitute;
import com.rtg.index.hash.ngs.instances.SubstituteIndel;
import com.rtg.util.test.RandomDna;
import com.rtg.util.Utils;
import com.rtg.util.integrity.Exam;

/**
 */
public class MaskTest extends AbstractSplitTest {

  @Override
  protected NgsHashFunction getHashFunction(final ReadCall readCall, final TemplateCall templateCall) {
    // not used
    return null;
  }


  protected void check(final Skeleton sk, final String dna, final String expected) throws IOException {
    final StringWriter sb = new StringWriter();
    final ReadCall rcall = new ReadCallMock(sb);
    final TemplateCall call = new TemplateCallMock(sb);
    final HashFunctionFactory f = Mask.factory(sk, false);
    final NgsHashFunction hf = f.create(rcall, call);
    //System.err.println("NgsHashFunction:" + hf);
    Exam.integrity(hf);
    sb.append(hf.toString()).append(LS).append("number windows=").append(String.valueOf(hf.numberWindows())).append(LS);
    hf.templateSet(1234, 23);
    encode(hf, dna);
    hf.readAll(5, false);
    hf.templateForward(7);
    assertEquals(expected, sb.toString());
  }

  protected void check(final Skeleton sk) throws IOException {
    final String dna = RandomDna.random(sk.readLength());
    final HashFunctionFactory f = Mask.factory(sk, false);
    final int indel = sk.indels();
    final Substitute sub = new Substitute(dna, f, true);
    sub.substituteProtected(sk.substitutions());
    if (indel > 0) {
      SubstituteIndel.checkIndel(f, dna, indel, 0/*cg*/, 7);
    }
  }

  public void testInvalid() {
    final Skeleton sk = new Skeleton(36, 36, 2, 0, 1);
    final String exp = ""
      + "Invalid skeleton:" + LS
      + "Inputs:invalid" + LS
      + "  Read length   36" + LS
      + "  Window length 36" + LS
      + "  Substitutions 2" + LS
      + "  Indels        0" + LS
      + "  Indel length  1" + LS
      ;
    try {
      Mask.factory(sk, false);
      fail();
    } catch (final Exception e) {
      assertEquals(exp, e.getMessage());
    }
  }
  public void test() throws IOException {
    final Skeleton sk = new Skeleton(36, 12, 2, 0, 1);
    check(sk);
  }

  private static final String SUFFIX = "done" + LS;

  private static final String EXPECTED_1 = ""
    + "Mask l=4 w=4 s=0 i=0" + LS
    + "number windows=1" + LS
    + "set name=1234 length=23" + LS
    + "readCall index=0 id=5" + LS
    + " hash=00000000:00000000:00000000:00000000:00000000:00000000:00000000:11001010" + LS
    + "templateCall position=7 index=0" + LS
    + " hash=00000000:00000000:00000000:00000000:00000000:00000000:00000000:11001010" + LS
    + SUFFIX
    ;

  public void test1() throws IOException {
    final Skeleton sk = new Skeleton(4, 4, 0, 0, 1);
    check(sk);
    final String dna = "TGCA";
    assertEquals(4, dna.length());
    check(sk, dna, EXPECTED_1);
  }

  private static final String EXPECTED_2 = ""
    + "Mask l=4 w=2 s=1 i=1" + LS
    + "number windows=2" + LS
    + "set name=1234 length=23" + LS
    + "readCall index=0 id=5" + LS
    + " hash=00000000:00000000:00000000:00000000:00000000:00000000:00000000:00000010" + LS
    + "readCall index=1 id=5" + LS
    + " hash=00000000:00000000:00000000:00000000:00000000:00000000:00000000:00001110" + LS
    + "templateCall position=7 index=0" + LS
    + " hash=00000000:00000000:00000000:00000000:00000000:00000000:00000000:00000010" + LS
    + "templateCall position=7 index=1" + LS
    + " hash=00000000:00000000:00000000:00000000:00000000:00000000:00000000:00001110" + LS
    + SUFFIX
    ;
  public void test2() throws IOException {
    final Skeleton sk = new Skeleton(4, 2, 1, 1, 1);
    check(sk);
    final String dna = "TGCA";
    assertEquals(4, dna.length());
    check(sk, dna, EXPECTED_2);
  }

  private static final String EXPECTED_3 = ""
    + "Mask l=36 w=18 s=2 i=1" + LS
    + "number windows=6" + LS
    + "set name=1234 length=23" + LS
    + "readCall index=0 id=5" + LS
    + " hash=00000000:00000000:00000000:00000000:00000000:00000011:11111110:00000000" + LS
    + "readCall index=1 id=5" + LS
    + " hash=00000000:00000000:00000000:00001111:11111000:00000000:00000000:00000000" + LS
    + "readCall index=2 id=5" + LS
    + " hash=00000000:00000000:00000000:00001111:11111000:00000011:11111110:00000000" + LS
    + "readCall index=3 id=5" + LS
    + " hash=00000000:00000000:00000000:00001111:11111000:00000000:00000001:11111111" + LS
    + "readCall index=4 id=5" + LS
    + " hash=00000000:00000000:00000000:00001111:11111000:00000011:11111111:11111111" + LS
    + "readCall index=5 id=5" + LS
    + " hash=00000000:00000000:00000000:00001111:11111111:11111111:11111110:00000000" + LS
    + "templateCall position=7 index=0" + LS
    + " hash=00000000:00000000:00000000:00000000:00000000:00000011:11111110:00000000" + LS
    + "templateCall position=7 index=1" + LS
    + " hash=00000000:00000000:00000000:00001111:11110000:00000000:00000010:00000000" + LS
    + "templateCall position=7 index=1" + LS
    + " hash=00000000:00000000:00000000:00001111:11111000:00000000:00000000:00000000" + LS
    + "templateCall position=7 index=1" + LS
    + " hash=00000000:00000000:00000000:00001111:11111000:00000010:00000000:00000000" + LS
    + "templateCall position=7 index=2" + LS
    + " hash=00000000:00000000:00000000:00001111:11111000:00000011:11111100:00000000" + LS
    + "templateCall position=7 index=2" + LS
    + " hash=00000000:00000000:00000000:00001111:11111000:00000011:11111110:00000000" + LS
    + "templateCall position=7 index=2" + LS
    + " hash=00000000:00000000:00000000:00000111:11111000:00000001:11111110:00000000" + LS
    + "templateCall position=7 index=3" + LS
    + " hash=00000000:00000000:00000000:00001111:11111000:00000000:00000001:11111111" + LS
    + "templateCall position=7 index=4" + LS
    + " hash=00000000:00000000:00000000:00001111:11111000:00000011:11111101:11111111" + LS
    + "templateCall position=7 index=4" + LS
    + " hash=00000000:00000000:00000000:00001111:11111000:00000011:11111111:11111111" + LS
    + "templateCall position=7 index=4" + LS
    + " hash=00000000:00000000:00000000:00000111:11111000:00000001:11111111:11111111" + LS
    + "templateCall position=7 index=5" + LS
    + " hash=00000000:00000000:00000000:00001111:11111111:11111111:11111110:00000000" + LS
    + "done" + LS
    ;
  public void test3() throws IOException {
    final Skeleton sk = new Skeleton(36, 18, 2, 1, 1);
    //System.err.println(sk.toString());
    check(sk);
    final String dna = "TTTTTTTTT" + "GGGGGGGGG" + "CCCCCCCCC" + "AAAAAAAAA";
    assertEquals(36, dna.length());
    check(sk, dna, EXPECTED_3);
  }

  public void testMaskFactory() {
    final Skeleton sk = new Skeleton(36, 18, 2, 1, 1);
    final HashFunctionFactory f = Mask.factory(sk, false);
    assertEquals(36, f.hashBits());
    assertEquals(6, f.numberWindows());
    assertEquals(36, f.windowBits());
  }

  public void testCGLAdjust() {
    assertEquals(0, Mask.cglAdjust(0L));
    assertEquals("00000000:00000000:00000000:00000111:11111111:11111111:11111111:11111111", Utils.toBitsSep(Mask.cglAdjust(Utils.fromBits("11111111:11111111:11111111:00000011:11111111"))));
    assertEquals("00000000:00000000:00000000:00000000:11100011:11111111:11111111:11111111", Utils.toBitsSep(Mask.cglAdjust(Utils.fromBits("00011000:11111111:11111111:00000011:11111111"))));
    assertEquals("00000000:00000000:00000000:00000000:01110011:11111111:11111111:11111111", Utils.toBitsSep(Mask.cglAdjust(Utils.fromBits("00001100:11111111:11111111:00000011:11111111"))));
    assertEquals("00000000:00000000:00000000:00000100:01100000:00000000:00000110:00000001", Utils.toBitsSep(Mask.cglAdjust(Utils.fromBits("10001000:00000000:00000001:00000010:00000001"))));

  }
}
