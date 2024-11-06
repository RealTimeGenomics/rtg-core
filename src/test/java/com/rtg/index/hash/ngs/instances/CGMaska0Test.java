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

import static com.rtg.util.StringUtils.LS;

import java.io.IOException;

import com.rtg.index.hash.ngs.NgsHashFunction;
import com.rtg.index.hash.ngs.ReadCall;
import com.rtg.index.hash.ngs.TemplateCall;

import org.junit.Assert;

/**
 */
public class CGMaska0Test extends AbstractSplitTest {

  @Override
  protected NgsHashFunction getHashFunction(final ReadCall readCall, final TemplateCall templateCall) {
    return new CGMaska0(readCall, templateCall) {
      @Override
      protected long hash(final long x) {
        return x;
      }
    };
  }

  private static final String PREFIX = ""
    + "CG l=35 w=18 e=0" + LS
    + "number windows=2" + LS
    + "set name=1234 length=23" + LS;

  private static final String SUFFIX = "done" + LS;

  private static final String EXPECTED_1 = ""
    + PREFIX
    + "readCall index=0 id=5" + LS
    + " hash=00000000:00000000:00000000:00001100:11001100:11001110:10101010:10101010" + LS
    + "readCall index=1 id=5" + LS
    + " hash=00000000:00000000:00000000:00000001:10011001:10011001:01010101:01010101" + LS
    + "templateCall position=3 index=0" + LS
    + " hash=00000000:00000000:00000000:00001100:11001100:11001110:10101010:10101010" + LS
    + "templateCall position=3 index=1" + LS
    + " hash=00000000:00000000:00000000:00000001:10011001:10011001:01010101:01010101" + LS
    + SUFFIX
    ;
  public void test1() throws IOException {
    final String dnar = "TGCATGCATGCATGCATGCATGCAT" +            "GCATGCATGC";
    final String dnat = "TGCATGCATGCATGCATGCATGCAT" + "aaaaaa" + "GCATGCATGC";
    assertEquals(35, dnar.length());
    assertEquals(41, dnat.length());
    check(dnar, dnat, EXPECTED_1);
  }

  private static final String EXPECTED_2 = ""
    + PREFIX
    + "readCall index=0 id=5" + LS
    + " hash=00000000:00000000:00000000:00001111:11111111:11111111:11111111:11111111" + LS
    + "readCall index=1 id=5" + LS
    + " hash=00000000:00000000:00000000:00000111:11111111:11111101:11111111:11111111" + LS
    + "templateCall position=3 index=0" + LS
    + " hash=00000000:00000000:00000000:00001111:11111111:11111111:11111111:11111111" + LS
    + "templateCall position=3 index=1" + LS
    + " hash=00000000:00000000:00000000:00000111:11111111:11111101:11111111:11111111" + LS
    + SUFFIX
    ;
  public void test2() throws IOException {
    final String dnar = "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT";
    final String dnat = "TTTTTTTTTTTTTTTTTTTTTTTTT" + "aaaaaa" + "TTTTTTTTTT";
    assertEquals(35, dnar.length());
    assertEquals(41, dnat.length());
    check(dnar, dnat, EXPECTED_2);
  }

  private static final String EXPECTED_3 = ""
    + PREFIX
    + "readCall index=0 id=5" + LS
    + " hash=00000000:00000000:00000000:00001111:11111111:11111100:00000000:00000000" + LS
    + "readCall index=1 id=5" + LS
    + " hash=00000000:00000000:00000000:00000111:11111111:11111100:00000000:00000000" + LS
    + "templateCall position=3 index=0" + LS
    + " hash=00000000:00000000:00000000:00001111:11111111:11111100:00000000:00000000" + LS
    + "templateCall position=3 index=1" + LS
    + " hash=00000000:00000000:00000000:00000111:11111111:11111100:00000000:00000000" + LS
    + SUFFIX
    ;
  public void test3() throws IOException {
    final String dnar = "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG";
    final String dnat = "GGGGGGGGGGGGGGGGGGGGGGGGG" + "aaaaaa" + "GGGGGGGGGG";
    assertEquals(35, dnar.length());
    assertEquals(41, dnat.length());
    check(dnar, dnat, EXPECTED_3);
  }

  public void testEmpty() throws IOException {
    check("", "", PREFIX);
  }

  public void testFactory() {
    assertEquals(36, CGMaska0.FACTORY.hashBits());
    assertEquals(36, CGMaska0.FACTORY.windowBits());
    assertEquals(2, CGMaska0.FACTORY.numberWindows());
    final NgsHashFunction hf = CGMaska0.FACTORY.create(new ReadCallMock(null), new TemplateCallMock(null));
    assertTrue(hf != null);
    assertTrue(hf instanceof CGMaska0);
    assertEquals(hf.numberWindows(), CGMaska0.FACTORY.numberWindows());
  }

  public void testAllSubstitutions() throws IOException {
    final String str = "acgacgtgacacccgtacgtaccccgtgacacccg";
    assertEquals(35, str.length());
    final SubstituteCG sub = new SubstituteCG(str, CGMaska0.FACTORY, true, 6, 0);
    Assert.assertEquals(0, sub.substituteProtected(1));
  }

  public void testIndel() throws IOException {
    final String str = "acgacgtgacacccgtacgtaccccgtgacacccg";
    assertEquals(35, str.length());
    SubstituteIndel.checkIndel(CGMaska0.FACTORY, str, 1, 6);
  }
}

