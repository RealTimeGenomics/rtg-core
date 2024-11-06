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

package com.rtg.variant.bayes;
import static com.rtg.util.StringUtils.LS;

import com.rtg.variant.bayes.snp.DescriptionCommon;
import com.rtg.variant.bayes.snp.DescriptionSnp;

import junit.framework.TestCase;

/**
 */
public class EvidenceTest extends TestCase {

  private static class MockEvidence extends Evidence {
    private final boolean mForward;
    MockEvidence(boolean forward) {
      super(DescriptionSnp.SINGLETON, 0.01);
      mForward = forward;
    }
    @Override
    public double error() {
      return 0.1;
    }

    @Override
    public int read() {
      return 0;
    }

    @Override
    public double probability(int index) {
      return (index + 1) / 10.0;
    }
    @Override
    public double pe() {
      return 0.2;
    }

    @Override
    public int getReadBasesLeft() {
      return 0;
    }

    @Override
    public int getReadBasesRight() {
      return 0;
    }

    @Override
    public void setReadBasesLeft(int readBasesLeft) {
      throw new UnsupportedOperationException();
    }

    @Override
    public void setReadBasesRight(int readBaseRight) {
      throw new UnsupportedOperationException();
    }

    @Override
    public boolean isReadPaired() {
      return false;
    }

    @Override
    public boolean isFirst() {
      return true;
    }

    @Override
    public boolean isMated() {
      return false;
    }

    @Override
    public boolean isUnmapped() {
      return false;
    }

    @Override
    public boolean isForward() {
      return mForward;
    }
  }

  public void test() {
    final Evidence di = new MockEvidence(true);
    di.globalIntegrity();
    assertTrue(di.description() instanceof DescriptionCommon);
    assertEquals(0.01, di.mapError());
    assertEquals(0, di.read());
    final String exp = ""
        + "Evidence > pe=0.200 error=0.100 notMap=0.010" + LS
        + " A -2.303" + LS
        + " C -1.609" + LS
        + " G -1.204" + LS
        + " T -0.916" + LS
        ;
    assertEquals(exp, di.toString());
  }

  public void testReverse() {
    final Evidence di = new MockEvidence(false);
    di.globalIntegrity();
    final String exp = ""
        + "Evidence < pe=0.200 error=0.100 notMap=0.010" + LS
        + " A -2.303" + LS
        + " C -1.609" + LS
        + " G -1.204" + LS
        + " T -0.916" + LS
        ;
    assertEquals(exp, di.toString());
  }
}
