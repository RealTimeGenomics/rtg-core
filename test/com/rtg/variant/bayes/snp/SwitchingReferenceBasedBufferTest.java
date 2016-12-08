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
