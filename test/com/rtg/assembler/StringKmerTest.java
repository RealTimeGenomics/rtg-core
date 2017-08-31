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

package com.rtg.assembler;


import com.rtg.assembler.graph.implementation.ContigString;
import com.rtg.util.diagnostic.Diagnostic;

/**
 */
public class StringKmerTest extends AbstractKmerTest {
  @Override
  protected StringKmer kmer(String str) {
    return new StringKmer(str);
  }

  public void testFactory() {
    Diagnostic.setLogStream();
    assertEquals("ACGTT", StringKmer.factory().make(new byte[] {1, 2, 3, 4, 4}, 0, 5).toString());
    assertEquals("GT", StringKmer.factory().make(new byte[] {1, 2, 3, 4, 4}, 2, 4).toString());
    assertEquals("GT", StringKmer.factory().make(new ContigString("ACGTT"), 2, 4).toString());
  }
}
