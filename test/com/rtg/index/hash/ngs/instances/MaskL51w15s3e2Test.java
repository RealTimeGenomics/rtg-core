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

/**
 */
public class MaskL51w15s3e2Test extends AbstractSplitTest {

  @Override
  protected NgsHashFunction getHashFunction(final ReadCall readCall, final TemplateCall templateCall) {
    return new MaskL51w15s3e2(readCall, templateCall) {
      @Override
      protected long hash(final long x) {
        return x;
      }
    };
  }

  public void testFactory() {
    assertEquals(30, MaskL51w15s3e2.FACTORY.hashBits());
    assertEquals(30, MaskL51w15s3e2.FACTORY.windowBits());
    assertEquals(9, MaskL51w15s3e2.FACTORY.numberWindows());
    final NgsHashFunction hf = MaskL51w15s3e2.FACTORY.create(new ReadCallMock(null), new TemplateCallMock(null));
    assertTrue(hf != null);
    assertTrue(hf instanceof MaskL51w15s3e2);
    assertEquals(hf.numberWindows(), MaskL51w15s3e2.FACTORY.numberWindows());
    assertEquals("MaskL51w15s3e2 l=51 w=15 s=3 e=2", hf.toString());
  }

  /**
   * Check that all 0, 1, 2 substitutions on the string are found.
   */
  public void testAllSubstitutions() throws IOException {
    final String str = "acgacgtgacacccgtacgtaccccgtgacacccgtccgtacgtaccccgt";
    assertEquals(51, str.length());
    final Substitute sub = new Substitute(str, MaskL51w15s3e2.FACTORY, true);
    sub.substituteProtected(3);

  }

  /**
   * Check that all 0, 1, 2 substitutions on the string are found.
   * @throws IOException
   */
  public void testIndel() throws IOException {
    final String str = "acgacgtgacacccgtacgtaccccgtgacacccgtccgtacgtaccccgt";
    assertEquals(51, str.length());
    SubstituteIndel.checkIndel(MaskL51w15s3e2.FACTORY, str, 2, 0/*cg*/, 7);
  }
}

