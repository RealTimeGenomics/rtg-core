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
package com.rtg.position.output;

import static com.rtg.mode.Protein.A;
import static com.rtg.mode.Protein.C;
import static com.rtg.mode.Protein.D;
import static com.rtg.mode.Protein.K;
import static com.rtg.mode.Protein.L;
import static com.rtg.mode.Protein.M;
import static com.rtg.mode.Protein.R;
import static com.rtg.mode.Protein.STOP;
import static com.rtg.mode.Protein.X;

import com.rtg.mode.Protein;
import com.rtg.util.StringUtils;

import junit.framework.TestCase;

/**
 */
public class Blosum62Test extends TestCase {

  /**
   * Test method for {@link Blosum62#score(int, int)}.
   */
  public final void testScore() {
    check(4, A, A);
    check(-1, A, R);
    check(1, STOP, STOP);
    check(-1, K, D);
    check(2, M, L);
    check(0, X, A);
    check(-1, X, R);
    check(-2, X, C);
    check(-4, X, STOP);
  }

  private void check(final int sc, final Protein a, final Protein b) {
    assertEquals(sc, Blosum62.SINGLETON.score(a.ordinal(), b.ordinal()));
    assertEquals(sc, Blosum62.SINGLETON.score(b.ordinal(), a.ordinal()));
  }

  public final void testToStringStringBuilder() {
    final String str = Blosum62.SINGLETON.toString();
    //System.err.println(str);
    assertTrue(str.startsWith("[-1, -4, 0, -1, -1, -1, -2,"));
    assertTrue(str.endsWith("-1, -2, -2, 0, -3, -1, 4]" + StringUtils.LS));
  }

}
