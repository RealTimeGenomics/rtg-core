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
package com.rtg.index.hash.ngs.protein;

import com.rtg.index.hash.ngs.ReadCall;
import com.rtg.index.hash.ngs.general.SingleMask;

import junit.framework.TestCase;

/**
 */
public class ProteinExtractReadTest extends TestCase {

  public void test() {
    final StringBuilder sb = new StringBuilder();
    final SingleMask ms = new SingleMask(4, 8, 0, 1);
    final ProteinExtractRead r = new ProteinExtractRead(ms, new ReadCall() {
        @Override
        public void readCall(final int id, final long hash, final int index) {
          if (hash == 0) {
            fail();
          }
        }
      }, 1);
    assertTrue(r.integrity());
    r.masked(0);
    r.toString(sb);
    final String s = sb.toString();
    assertTrue(s.contains(" index=1 id=-1"));
    assertTrue(s.contains("ProteinExtractRead readCall="));
  }

}
