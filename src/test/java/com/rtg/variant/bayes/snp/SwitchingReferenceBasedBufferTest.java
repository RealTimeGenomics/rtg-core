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
package com.rtg.variant.bayes.snp;

import com.rtg.variant.bayes.ReferenceBasedFactory;

import junit.framework.TestCase;

/**
 */
public class SwitchingReferenceBasedBufferTest extends TestCase {

  private static class MockModel {
    private final String mName;

    MockModel(String prefix, final int id, final int ref) {
      mName = prefix + ":" + id + ":" + ref;
    }

    @Override
    public String toString() {
      return mName;
    }
  }

  private static class Fac implements ReferenceBasedFactory<MockModel> {
    private int mCount = 0;
    private final String mPrefix;
    Fac(String prefix) {
      mPrefix = prefix;
    }
    @Override
    public MockModel make(final int ref) {
      return new MockModel(mPrefix, mCount++, ref);
    }
  }

  public void test() {
    final byte[] template = {1, 2, 3, 4, 3, 2, 1, 2, 0, 4, 3, 2, 1, 2, 3, 4, 3, 2, 1, 2, 3, 4, 3, 2, 1, 2, 3, 4, 3, 2, 1};
    final ReferenceBasedBuffer<MockModel> cb = new SwitchingReferenceBasedBuffer<>(1, new Fac("A"), new Fac("B"), 9, template, 0);

    assertEquals(0, cb.base());
    assertEquals("A:0:3", cb.get(3).toString());
    assertEquals(0, cb.base());

    final MockModel step = cb.step();
    assertEquals(1, cb.base());
    assertEquals("A:1:0", step.toString());

    assertEquals(1, cb.find(1));
    assertEquals(9, cb.find(9));
    assertEquals(0, cb.find(10));

    assertEquals("B:0:3", cb.get(9).toString()); // First object made after the switch
    assertEquals("A:2:1", cb.get(1).toString());
    assertEquals("A:3:-1", cb.get(8).toString());

    assertEquals(10, cb.find(11));

    for (int i = 1; i < 24; ++i) {
      final MockModel model = cb.step();
      assertNotNull("" + i, model);
      cb.globalIntegrity();
    }
    assertEquals(23, cb.find(24));
    assertEquals(0, cb.find(25));
    cb.step();
    badFind(cb, 1, 25);
  }

  private void badFind(final ReferenceBasedBuffer<?> mb, final int f, final int b) {
    try {
      mb.find(f);
      fail();
    } catch (final RuntimeException e) {
      assertEquals("Index less than base. index=" + Integer.toString(f) + " base=" + Integer.toString(b), e.getMessage());
    }
  }
}
