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

import com.rtg.index.hash.ngs.NgsHashFunction;
import com.rtg.index.hash.ngs.ReadCall;
import com.rtg.index.hash.ngs.TemplateCall;

import org.junit.Assert;

/**
 */
public class CG2Maska1Test extends AbstractSplitTest {

  @Override
  protected NgsHashFunction getHashFunction(final ReadCall readCall, final TemplateCall templateCall) {
    return new CG2Maska1(readCall, templateCall) {
      @Override
      protected long hash(final long x) {
        return x;
      }
    };
  }

  public void testFactory() {
    assertEquals(30, CG2Maska1.FACTORY.hashBits());
    assertEquals(30, CG2Maska1.FACTORY.windowBits());
    assertEquals(2, CG2Maska1.FACTORY.numberWindows());
    final NgsHashFunction hf = CG2Maska1.FACTORY.create(new ReadCallMock(null), new TemplateCallMock(null));
    assertTrue(hf != null);
    assertTrue(hf instanceof CG2Maska1);
    assertEquals(hf.numberWindows(), CG2Maska1.FACTORY.numberWindows());
  }

  public void testOverlap4() throws IOException {
    final String dnar = "TGCATGCATG"
                            + "CATGCATGCATGCATGCAT";
    final String dnat = "TGCATGCATGCATGCATGCATGCAT";
    assertEquals(29, dnar.length());
    assertEquals(25, dnat.length());
    checkNano("cg2maska1-overlap4.txt", dnar, dnat);
  }

  public void testOverlap3() throws IOException {
    final String dnar = "TGCATGCATG"
                             + "ATGCATGCATGCATGCATG";
    final String dnat = "TGCATGCATGCATGCATGCATGCATG";
    assertEquals(29, dnar.length());
    assertEquals(26, dnat.length());
    checkNano("cg2maska1-overlap3.txt", dnar, dnat);
  }

  public void test2() throws IOException {
    final String dnar = "TTTTTTTTTTTTTTTTTTTTTTTTTTTTT";
    final String dnat = "TTTTTTTTTTTTTTTTTTTTTTTTT";
    assertEquals(29, dnar.length());
    assertEquals(25, dnat.length());
    checkNano("cg2maska1-test2.txt", dnar, dnat);
  }

  public void test3() throws IOException {
    final String dnar = "GGGGGGGGGGGGGGGGGGGGGGGGGGGGG";
    final String dnat = "GGGGGGGGGGGGGGGGGGGGGGGGG";
    assertEquals(29, dnar.length());
    assertEquals(25, dnat.length());
    checkNano("cg2maska1-test3.txt", dnar, dnat);
  }

  public void testEmpty() throws IOException {
    checkNano("cg2maska1-empty.txt", "", "");
  }


  public void testAllSubstitutions() throws IOException {
    final String str = "tgacacccgt"
                           + "ccgtaccccgtgacacccg";
    assertEquals(29, str.length());
    final SubstituteCG sub = new SubstituteCG(str, CG2Maska1.FACTORY, true, 0, 4);
    final int missed1 = sub.substituteProtected(1);
    Assert.assertEquals(sub.toString(), 4, missed1);
    final int missed2 = sub.substituteProtected(2);
    Assert.assertEquals(sub.toString(), 220,  missed2);
  }

//  public void testIndel() throws IOException {
//    final String str = "acgacgtgacacccgtacgtaccccgtgacacccg";
//    assertEquals(29, str.length());
//    SubstituteIndel.checkIndel(CG2Maska0.FACTORY, str, 1);
//  }
}

