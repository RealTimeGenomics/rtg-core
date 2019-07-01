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
