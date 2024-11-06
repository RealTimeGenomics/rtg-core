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
