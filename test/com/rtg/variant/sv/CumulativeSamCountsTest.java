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

package com.rtg.variant.sv;

import static com.rtg.util.StringUtils.LS;

import java.io.ByteArrayInputStream;
import java.io.IOException;


/**
 */
public class CumulativeSamCountsTest extends AbstractSamCountsTest {

  @Override
  protected SamCounts getCounts(int length) {
    return new CumulativeSamCounts(length, null);
  }

  public void testCorrections() throws IOException {
    final String coStr = ""
      + "seq1 0 1.0" + LS
      + "seq1 1 2.0" + LS
      + "seq1 2 0.25" + LS
      + "seq1 5 1.5" + LS
      ;
    final Corrections co = new Corrections(new ByteArrayInputStream(coStr.replace(' ', '\t').getBytes()), 1);
    final CumulativeSamCounts csc = new CumulativeSamCounts(10, co);
    csc.increment(0);
    csc.increment(0);
    csc.increment(1);
    csc.increment(3);
    assertEquals(2.0, csc.count(0, 0));
    assertEquals(0.5, csc.count(0, 1));
    assertEquals(0.0, csc.count(0, 3));
  }
}
