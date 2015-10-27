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

import com.rtg.index.hash.ngs.NgsHashFunction;
import com.rtg.index.hash.ngs.ReadCall;
import com.rtg.index.hash.ngs.TemplateCall;

import junit.framework.Assert;

/**
 */
public class CG2Maska15Test extends AbstractSplitTest {

  @Override
  protected NgsHashFunction getHashFunction(final ReadCall readCall, final TemplateCall templateCall) {
    return new CG2Maska15(readCall, templateCall) {
      @Override
      protected long hash(final long x) {
        return x;
      }
    };
  }

  public void testFactory() {
    assertEquals(30, CG2Maska15.FACTORY.hashBits());
    assertEquals(30, CG2Maska15.FACTORY.windowBits());
    assertEquals(6, CG2Maska15.FACTORY.numberWindows());
    final NgsHashFunction hf = CG2Maska15.FACTORY.create(new ReadCallMock(null), new TemplateCallMock(null));
    assertTrue(hf != null);
    assertTrue(hf instanceof CG2Maska15);
    assertEquals(hf.numberWindows(), CG2Maska15.FACTORY.numberWindows());
  }

  public void testOverlap3() throws IOException {
    final String dnar = "TGCATGCATG"
      + "ATGCATGCATGCATGCATG";
    final String dnat = "TGCATGCATGCATGCATGCATGCATG";
    assertEquals(29, dnar.length());
    assertEquals(26, dnat.length());
    checkNano("cg2maska15-overlap3.txt", dnar, dnat);
  }

  public void testOverlap4() throws IOException {
    final String dnar = "TGCATGCATG"
                            + "CATGCATGCATGCATGCAT";
    final String dnat = "TGCATGCATGCATGCATGCATGCAT";
    assertEquals(29, dnar.length());
    assertEquals(25, dnat.length());
    checkNano("cg2maska15-overlap4.txt", dnar, dnat);
  }

  public void testOverlap5() throws IOException {
    final String dnar = "TGCATGCATG"
      + "GCATGCATGCATGCATGCA";
    final String dnat = "TGCATGCATGCATGCATGCATGCA";
    assertEquals(29, dnar.length());
    assertEquals(24, dnat.length());
    checkNano("cg2maska15-overlap5.txt", dnar, dnat);
  }

  public void test2() throws IOException {
    final String dnar = "TTTTTTTTTTTTTTTTTTTTTTTTTTTTT";
    final String dnat = "TTTTTTTTTTTTTTTTTTTTTTTTT";
    assertEquals(29, dnar.length());
    assertEquals(25, dnat.length());
    checkNano("cg2maska15-test2.txt", dnar, dnat);
  }

  public void test3() throws IOException {
    final String dnar = "GGGGGGGGGGGGGGGGGGGGGGGGGGGGG";
    final String dnat = "GGGGGGGGGGGGGGGGGGGGGGGGG";
    assertEquals(29, dnar.length());
    assertEquals(25, dnat.length());
    checkNano("cg2maska15-test3.txt", dnar, dnat);
  }

  public void testEmpty() throws IOException {
    checkNano("cg2maska15-empty.txt", "", "");
  }


  public void testAllSubstitutions() throws IOException {
    final String str = "tgacacccgt"
                           + "ccgtaccccgtgacacccg";
    assertEquals(29, str.length());
    final SubstituteCG sub = new SubstituteCG(str, CG2Maska15.FACTORY, true, 0, 4);
    final int missed1 = sub.substituteProtected(1);
    Assert.assertEquals(sub.toString(), 0, missed1);
    final int missed2 = sub.substituteProtected(2);
    Assert.assertEquals(sub.toString(), 94, missed2);
  }

//  public void testIndel() throws IOException {
//    final String str = "acgacgtgacacccgtacgtaccccgtgacacccg";
//    assertEquals(29, str.length());
//    SubstituteIndel.checkIndel(CG2Maska0.FACTORY, str, 1);
//  }
}

