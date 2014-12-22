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
package com.rtg.mode;

import java.util.Arrays;

import com.rtg.util.Utils;

import junit.framework.TestCase;

/**
 * Tests the corresponding class.
 *
 */
public class DnaUtilsTest extends TestCase {

  public void testToCodes01() {
    assertEquals("AAA", DnaUtils.toCodes(Utils.fromBits(""), Utils.fromBits(""), 3));
    assertEquals("TGC", DnaUtils.toCodes(Utils.fromBits("110"), Utils.fromBits("101"), 3));
    assertEquals("AT", DnaUtils.toCodes(Utils.fromBits("101"), Utils.fromBits("101"), 2));
    assertEquals("T", DnaUtils.toCodes(Utils.fromBits("101"), Utils.fromBits("101"), 1));
    assertEquals("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATAT", DnaUtils.toCodes(Utils.fromBits("101"), Utils.fromBits("101"), 64));
    try {
      DnaUtils.toCodes(0, 0);
    } catch (final RuntimeException e) {
      assertEquals("length out of range=0", e.getMessage());
    }
    try {
      DnaUtils.toCodes(0, 65);
    } catch (final RuntimeException e) {
      assertEquals("length out of range=65", e.getMessage());
    }
  }

  public void testToCodes() {
    assertEquals("AAA", DnaUtils.toCodes(Utils.fromBits(""), 3));
    assertEquals("CAC", DnaUtils.toCodes(Utils.fromBits("101"), 3));
    assertEquals("AT", DnaUtils.toCodes(Utils.fromBits("101"), 2));
    assertEquals("G", DnaUtils.toCodes(Utils.fromBits("110"), 1));
    assertEquals("AAAAAAAAAAAAAAAAAAAAAAAAAAAAACAC", DnaUtils.toCodes(Utils.fromBits("101"), 32));
    try {
      DnaUtils.toCodes(0, 0);
    } catch (final RuntimeException e) {
      assertEquals("length out of range=0", e.getMessage());
    }
    try {
      DnaUtils.toCodes(0, 33);
    } catch (final RuntimeException e) {
      assertEquals("length out of range=33", e.getMessage());
    }
  }

  public void testReverseComplement() {
    assertEquals("", DnaUtils.reverseComplement(""));
    assertEquals("ncagt", DnaUtils.reverseComplement("actgn"));
    assertEquals("NCAGT", DnaUtils.reverseComplement("ACTGN"));
    assertEquals("tgcatgca", DnaUtils.reverseComplement("tgcatgca"));
  }

  public void testPhredToFastq() {
    final byte[] qualityScores = new byte[35];
    for (int i = 0; i < qualityScores.length; i++) {
      qualityScores[i] = (byte) i;
    }
    assertEquals("" + ((char) (4 + 33)), DnaUtils.phredToFastq(qualityScores, 4, 1, false));
    assertEquals("" + ((char) (31 + 33)) + ((char) (30 + 33)), DnaUtils.phredToFastq(qualityScores, 30, 2, true));
  }

  public void testEncodeString() {
    final String exp = "NACGT";
    final byte[] b = DnaUtils.encodeString(exp);
    final String s = DnaUtils.bytesToSequenceIncCG(b);
    assertEquals(exp, s);
  }

  public void testByteRevCo() {
    final byte[] src = {'A', 'C', 'C', 'T', 'G', 'G', 'A'};
    final byte[] dest = new byte[src.length];

    DnaUtils.reverseComplement(src, dest, src.length);
    assertTrue(Arrays.toString(dest), Arrays.equals(new byte[] {'T', 'C', 'C', 'A', 'G', 'G', 'T'}, dest));

    final byte[] dest2 = new byte[src.length];

    DnaUtils.reverseComplement(src, dest2, 5);
    assertTrue(Arrays.toString(dest2), Arrays.equals(new byte[] {'C', 'A', 'G', 'G', 'T', 0, 0}, dest2)); //cheating by using N as 0, oh well.
  }

}

