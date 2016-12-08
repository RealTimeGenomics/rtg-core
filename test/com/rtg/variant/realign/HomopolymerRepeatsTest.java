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
package com.rtg.variant.realign;

import java.util.Arrays;

import com.rtg.mode.DnaUtils;

import junit.framework.TestCase;

/**
 * Test class
 */
public class HomopolymerRepeatsTest extends TestCase {

  private static class ByteArray extends ByteArrayAdaptor {
    private final byte[] mBytes;

    /**
     * @param bytes underlying byte array being proxied.
     */
    ByteArray(byte[] bytes) {
      mBytes = bytes;
    }

    @Override
    public byte get(int index) {
      return mBytes[index];
    }

    @Override
    public int length() {
      return mBytes.length;
    }

  }

  private static final String TEMPLATE = "CAAATGGCA";
  private static final int[] EXP_FORWARD = {1, 0, 0, 3, 1, 0, 2, 1, 1};
  private static final int[] EXP_REVERSE = {1, 3, 0, 0, 1, 2, 0, 1, 1};

  //simple example
  public void test() {
    final HomopolymerRepeats counts = new HomopolymerRepeats(new ByteArray(DnaUtils.encodeString(TEMPLATE)));
    check(counts, EXP_FORWARD, EXP_REVERSE);
  }

  //all 1 long
  public void test2() {
    final String template = "ACGTGCAGTCGATCGACGTGTGTGTCAC";
    final HomopolymerRepeats counts = new HomopolymerRepeats(new ByteArray(DnaUtils.encodeString(template)));
    final int[] exp = new int[template.length()];
    Arrays.fill(exp, 1);
    check(counts, exp, exp);
  }

  private static final int[] EXP_FORWARD1 = {0, 2, 1, 0, 2, 1};
  private static final int[] EXP_REVERSE1 = {2, 0, 1, 2, 0, 1};

  //non-zero start and end
  public void test1() {
    final HomopolymerRepeats counts = new HomopolymerRepeats(new ByteArray(DnaUtils.encodeString(TEMPLATE)), 2, 8);
    check(counts, EXP_FORWARD1, EXP_REVERSE1);
  }


  private void check(final HomopolymerRepeats counts, final int[] forward, final int[] reverse) {
    assertEquals(forward.length, reverse.length);
    for (int i = 0; i < forward.length; ++i) {
      assertEquals(counts.forward(i), forward[i]);
      assertEquals(counts.reverse(i), reverse[i]);
    }
  }
}
