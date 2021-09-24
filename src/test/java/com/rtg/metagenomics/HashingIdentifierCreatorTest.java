/*
 * Copyright (c) 2017. Real Time Genomics Limited.
 *
 * Use of this source code is bound by the Real Time Genomics Limited Software Licence Agreement
 * for Academic Non-commercial Research Purposes only.
 *
 * If you did not receive a license accompanying this file, a copy must first be obtained by email
 * from support@realtimegenomics.com.  On downloading, using and/or continuing to use this source
 * code you accept the terms of that license agreement and any amendments to those terms that may
 * be made from time to time by Real Time Genomics Limited.
 */

package com.rtg.metagenomics;

import com.rtg.mode.DnaUtils;

import junit.framework.TestCase;

/**
 * Tests the corresponding class.
 */
public class HashingIdentifierCreatorTest extends TestCase {

  public void test() {
    final String read = "ACCCCG";
    final byte[] bases = read.getBytes();
    final byte[] rc = DnaUtils.reverseComplement(read).getBytes();
    assertEquals(HashingIdentifierCreator.irvineHash(bases), HashingIdentifierCreator.irvineHashRC(rc));
    final HashingIdentifierCreator creator = new HashingIdentifierCreator();
    assertEquals(-2373396206520789191L, creator.getIdentifier("read", bases, false));
    assertEquals(-2373396206520789191L, creator.getIdentifier("read", rc, true));
  }
}
