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

package com.rtg.index;

import com.rtg.AbstractTest;

/**
 * Test
 */
public class ProportionalRepeatFrequencyFilterMethodTest extends AbstractTest {

  public void testProportional() {
    final ProportionalRepeatFrequencyFilterMethod m = new ProportionalRepeatFrequencyFilterMethod(50, 100, 0);
    //discard half the hashes
    final SparseFrequencyHistogram sph = new SparseFrequencyHistogram();
    for (int i = 1; i <= 10; ++i) {
      sph.add(i, 10 / i);
    }
    m.internalInitializeProportional(sph, 100);
    assertTrue(m.keepHash(0, 2));
    assertTrue(m.keepHash(0, 3));
    assertTrue(m.keepHash(0, 4));
    assertFalse(m.keepHash(0, 5));
    assertFalse(m.keepHash(0, 6));
    assertFalse(m.keepHash(0, 9));
  }
}
