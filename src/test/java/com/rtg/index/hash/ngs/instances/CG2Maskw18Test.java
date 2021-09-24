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

import org.junit.Assert;

/**
 */
public class CG2Maskw18Test extends AbstractSplitTest {

  @Override
  protected NgsHashFunction getHashFunction(final ReadCall readCall, final TemplateCall templateCall) {
    return new CG2Maskw18(readCall, templateCall) {
      @Override
      protected long hash(final long x) {
        return x;
      }
    };
  }

  public void testFactory() {
    assertEquals(36, CG2Maskw18.FACTORY.hashBits());
    assertEquals(36, CG2Maskw18.FACTORY.windowBits());
    assertEquals(3, CG2Maskw18.FACTORY.numberWindows());
    final NgsHashFunction hf = CG2Maskw18.FACTORY.create(new ReadCallMock(null), new TemplateCallMock(null));
    assertTrue(hf != null);
    assertTrue(hf instanceof CG2Maskw18);
    assertEquals(hf.numberWindows(), CG2Maskw18.FACTORY.numberWindows());
  }

  public void testOverlap4() throws IOException {
    final String dnar = "TGCATGCATG"
                            + "CATGCATGCATGCATGCAT";
    final String dnat = "TGCATGCATGCATGCATGCATGCAT";
    assertEquals(29, dnar.length());
    assertEquals(25, dnat.length());
    checkNano("cg2maskw18-overlap4.txt", dnar, dnat);
  }

  public void testOverlap3() throws IOException {
    final String dnar = "TGCATGCATG"
                             + "ATGCATGCATGCATGCATG";
    final String dnat = "TGCATGCATGCATGCATGCATGCATG";
    assertEquals(29, dnar.length());
    assertEquals(26, dnat.length());
    checkNano("cg2maskw18-overlap3.txt", dnar, dnat);
  }

  public void test2() throws IOException {
    final String dnar = "TTTTTTTTTTTTTTTTTTTTTTTTTTTTT";
    final String dnat = "TTTTTTTTTTTTTTTTTTTTTTTTT";
    assertEquals(29, dnar.length());
    assertEquals(25, dnat.length());
    checkNano("cg2maskw18-test2.txt", dnar, dnat);
  }

  public void test3() throws IOException {
    final String dnar = "GGGGGGGGGGGGGGGGGGGGGGGGGGGGG";
    final String dnat = "GGGGGGGGGGGGGGGGGGGGGGGGG";
    assertEquals(29, dnar.length());
    assertEquals(25, dnat.length());
    checkNano("cg2maskw18-test3.txt", dnar, dnat);
  }

  public void testEmpty() throws IOException {
    checkNano("cg2maskw18-empty.txt", "", "");
  }


  public void testAllSubstitutions() throws IOException {
    final String str = "tgacacccgt"
                           + "ccgtaccccgtgacacccg";
    assertEquals(29, str.length());
    final SubstituteCG sub = new SubstituteCG(str, CG2Maskw18.FACTORY, true, 0, 4);
    final int missed1 = sub.substituteProtected(1);
    //System.err.println(sub.details());
    Assert.assertEquals(sub.toString(), 4, missed1);
    final int missed2 = sub.substituteProtected(2);
    //System.err.println(sub.details());
    Assert.assertEquals(sub.toString(), 257,  missed2);
  }

//  public void testIndel() throws IOException {
//    final String str = "acgacgtgacacccgtacgtaccccgtgacacccg";
//    assertEquals(29, str.length());
//    SubstituteIndel.checkIndel(CG2Maska0.FACTORY, str, 1);
//  }
}

