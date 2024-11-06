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

package com.rtg.variant.sv;


/**
 */
public class SamArrayTest extends AbstractSamCountsTest {

  @Override
  protected SamCounts getCounts(int length) {
    return new SamArray(length);
  }

  @Override
  protected void check10Inc(final SamCounts sa) {
    sa.increment(0);
    ((SamArray) sa).increment(4, 2.0);
  }

  public void testReverse() {
    final SamArray sa = new SamArray(5);
    for (int i = 0; i < 5; ++i) {
      sa.increment(i, i + 1);
    }
    assertEquals("[" + 5.0 + ", " + 5.0 + ", " + 4.0 + ", " + 3.0 + ", " + 2.0 + "]", sa.reverse(3, 0).toString());
    assertEquals("[" + 5.0 + ", " + 4.0 + ", " + 3.0 + ", " + 2.0 + ", " + 1.0 + "]", sa.reverse(3, 1).toString());
    assertEquals("[" + 4.0 + ", " + 3.0 + ", " + 2.0 + ", " + 1.0 + ", " + 1.0 + "]", sa.reverse(3, 2).toString());
  }
}
