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

package com.rtg.variant.bayes.multisample.family;

import com.rtg.util.MathUtils;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;
import com.rtg.variant.bayes.Code;

/**
 * Implements the parents haploid / diploid and children either haploid or diploid.
 */
public abstract class AbstractSegregationHaploid extends IntegralAbstract implements ISegregationScore {

  private final int[] mCounts = new int[2];
  protected final Code mCode;
  protected final UniqueId mUid;

  private static final class XYX extends AbstractSegregationHaploid {
    private XYX(Code code, int father, int mother) {
      super(code, father, mother);
    }
    @Override
    protected int childIndex(int child, boolean diploid) {
      final int ch0 = mUid.id(mCode.a(child));
      final int ch1 = mUid.id(mCode.bc(child));
      if (ch0 < 0 || ch1 < 0) {
        return -1;
      }
      if (ch0 == 0 && ch1 == 0) {
        return 0;
      } else if (ch0 == 0 && ch1 == 1) {
        return 1;
      } else if (ch0 == 1 && ch1 == 1) {
        if (diploid) {
          return -1;
        }
        return 1;
      }
      throw new RuntimeException();
    }

    @Override
    public boolean integrity() {
      super.integrity();
      Exam.assertEquals(2, mUid.numberIdsSoFar());
      return true;
    }
  }

  private static final class XYY extends AbstractSegregationHaploid {
    private XYY(Code code, int father, int mother) {
      super(code, father, mother);
    }
    @Override
    protected int childIndex(int child, boolean diploid) {
      final int ch0 = mUid.id(mCode.a(child));
      final int ch1 = mUid.id(mCode.bc(child));
      if (ch0 < 0 || ch1 < 0) {
        return -1;
      }
      if (ch0 == 0 && ch1 == 0) {
        if (diploid) {
          return -1;
        }
        return 0;
      } else if (ch0 == 0 && ch1 == 1) {
        return 0;
      } else if (ch0 == 1 && ch1 == 1) {
        return 1;
      }
      throw new RuntimeException();
    }

    @Override
    public boolean integrity() {
      super.integrity();
      Exam.assertEquals(2, mUid.numberIdsSoFar());
      return true;
    }
  }

  private static final class XYZ extends AbstractSegregationHaploid {
    private XYZ(Code code, int father, int mother) {
      super(code, father, mother);
    }
    @Override
    protected int childIndex(int child, boolean diploid) {
      final int ch0 = mUid.id(mCode.a(child));
      final int ch1 = mUid.id(mCode.bc(child));
      if (ch0 < 0 || ch1 < 0) {
        return -1;
      }
      if (ch0 == 0 && ch1 == 1) {
        return -1;
      }
      if (ch0 == 0 || ch1 == 0) {
        return 0;
      } else if (ch0 == 1 || ch1 == 1) {
        return 1;
      }
      return -1;
    }

    @Override
    public boolean integrity() {
      super.integrity();
      Exam.assertEquals(3, mUid.numberIdsSoFar());
      return true;
    }
  }

  /**
   * Get a segregation score instance that expects the father to be haploid.
   * @param code the code to use (should have range size of the number of alternative alleles plus one for the reference)
   * @param father the father hypothesis code (must be haploid)
   * @param mother the mother hypothesis code
   * @return segregation score for haploid parent implementation
   */
  public static ISegregationScore getHaploidInstance(Code code, int father, int mother) {
    assert code.homozygous(father);
    if (code.homozygous(mother)) {
      return new SegregationTrivial();
    }
    final UniqueId uid = new UniqueId(code.rangeSize());
    final int mo0 = uid.addId(code.a(mother));
    assert mo0 == 0;
    final int mo1 = uid.addId(code.bc(mother));
    assert mo1 == 1;
    final int fa0 = uid.addId(code.a(father));
    if (fa0 == 0) {
      return new XYX(code, father, mother);
    }
    if (fa0 == 1) {
      return new XYY(code, father, mother);
    }
    if (fa0 == 2) {
      return new XYZ(code, father, mother);
    }
    throw new RuntimeException();
  }

  private AbstractSegregationHaploid(Code code, int father, int mother) {
    mCode = code;
    mUid = new UniqueId(mCode.rangeSize());
    final int mo0 = mUid.addId(mCode.a(mother));
    final int mo1 = mUid.addId(mCode.bc(mother));
    assert mo0 == 0 && mo1 == 1;
    assert mCode.homozygous(father);
    final int fa1 = mUid.addId(mCode.a(father));
    assert fa1 <= 2;
  }

  @Override
  public void increment(int child, boolean diploid) {
    final int index = childIndex(child, diploid);
    if (index >= 0) {
      mCounts[index]++;
    }
  }

  @Override
  public double lnProbability() {
    final int total = mCounts[0] + mCounts[1];
    final double ret = MathUtils.logFactorial(total) - MathUtils.logFactorial(mCounts[0]) - MathUtils.logFactorial(mCounts[1]) - MathUtils.LOG_2 * total;
    assert !Double.isNaN(ret);
    return ret;
  }

  protected abstract int childIndex(int child, boolean diploid);

  @Override
  public boolean integrity() {
    Exam.assertTrue(mCounts[0] >= 0);
    Exam.assertTrue(mCounts[1] >= 0);
    Exam.assertTrue(mCounts.length == 2);
    Exam.assertTrue(mUid.numberIdsSoFar() <= 3);
    return true;
  }
}
