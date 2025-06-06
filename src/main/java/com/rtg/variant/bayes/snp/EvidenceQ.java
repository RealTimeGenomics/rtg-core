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

import com.rtg.util.integrity.Exam;
import com.rtg.variant.bayes.Code;
import com.rtg.variant.bayes.CodeDiploid;
import com.rtg.variant.bayes.Description;
import com.rtg.variant.bayes.Evidence;
import com.rtg.variant.util.arithmetic.LogApproximatePossibility;
import com.rtg.variant.util.arithmetic.PossibilityArithmetic;

/**
 * Generates a distribution from a primary prediction and a quality value.
 */
public final class EvidenceQ extends Evidence {

  private final int mRead;

  private final ThreadLocal<Integer> mReadBasesLeft = new ThreadLocal<>();

  private final ThreadLocal<Integer> mReadBasesRight = new ThreadLocal<>();

  private final double mQ;

  private final double mQt;

  private final double mQm;

  private final boolean mForward;
  private final boolean mReadPaired;
  private final boolean mFirst;
  private final boolean mMated;
  private final boolean mUnmapped;

  private final double[] mEvidentialPossibility;

  private static final PossibilityArithmetic ARITHMETIC = LogApproximatePossibility.SINGLETON;

  /**
   * @param description of the underlying haploid hypotheses.
   * @param read the hypothesis corresponding to the read.
   * @param r the probability that the read is not mapped to this position.
   * @param q the total probability that the hypothesis is not equal to <code>read</code>.
   * @param isForward true iff mapped in forward frame
   * @param isReadPaired true iff mapping is from a paired read
   * @param isFirst true iff the read if unpaired or the first of a pair
   * @param isMated true iff mapping from a paired read which is mated
   * @param isUnmapped if mapping from an unmapped arm
   */
  public EvidenceQ(final Description description, final int read, final double r, final double q, boolean isForward, boolean isReadPaired, boolean isFirst, boolean isMated, boolean isUnmapped) {
    super(description, r);
    mQ = q;
    mQt = q / (description.size() - 1);
    mQm = 1.0 - q;
    mRead = read;
    mForward = isForward;
    mReadPaired = isReadPaired;
    mFirst = isFirst;
    mMated = isMated;
    mUnmapped = isUnmapped;
    final Code code = new CodeDiploid(description().size());
    mEvidentialPossibility = new double[code.size()];
    final double rc = 1.0 - r;
    final double pE = pe();
    assert pE >= 0 && pE <= 1;
    final double pEr = r * pE;
    for (int i = 0; i < code.size(); ++i) {
      final double prob = 0.5 * (probability(code.a(i)) + probability(code.bc(i)));
      if (prob <= 0.0) {
        mEvidentialPossibility[i] = 1;
      } else {
        final double pr = prob * rc + pEr;
        mEvidentialPossibility[i] = ARITHMETIC.prob2Poss(pr);
      }
    }
  }
  /**
   * @param description of the underlying haploid hypotheses.
   * @param read the hypothesis corresponding to the read.
   * @param readBasesLeft bases to left of <code>read</code> on read
   * @param readBasesRight bases to right of <code>read</code> on read
   * @param r the probability that the read is not mapped to this position.
   * @param q the total probability that the reference is not equal to <code>read</code>.
   * @param isForward true if mapped in forward frame, false otherwise
   * @param isReadPaired true if mapping is from a paired read, false otherwise
   * @param isFirst true if mapping is from the first arm of a paired read
   * @param isMated true if mapping from a paired read which is mated, false otherwise
   * @param isUnmapped true if mapping from an unmapped arm
   */
  public EvidenceQ(final Description description, final int read, int readBasesLeft, int readBasesRight, final double r, final double q, boolean isForward, boolean isReadPaired, boolean isFirst, boolean isMated, boolean isUnmapped) {
    this(description, read, r, q, isForward, isReadPaired, isFirst, isMated, isUnmapped);
    mReadBasesLeft.set(readBasesLeft);
    mReadBasesRight.set(readBasesRight);
  }

  @Override
  public double probability(int index) {
    return index == mRead ? mQm : mQt;
  }

  @Override
  public double pe() {
    return 0.25;
  }

  /**
   * Precomputed probability to increment model with for this evidence.
   * @param hyp hypothesis
   * @return probability
   */
  public double logEvidentialProbability(final int hyp) {
    return mEvidentialPossibility[hyp];
  }

  @Override
  public double error() {
    return mQ;
  }

  @Override
  public int read() {
    return mRead;
  }

  /**
   * Be very careful when using this field, it is stored in a {@code ThreadLocal} so we can use the singleton cache
   */
  @Override
  public int getReadBasesLeft() {
    return mReadBasesLeft.get();
  }

  /**
   * Be very careful when using this field, it is stored in a {@code ThreadLocal} so we can use the singleton cache
   */
  @Override
  public int getReadBasesRight() {
    return mReadBasesRight.get();
  }

  /**
   * Be very careful when using this field, it is stored in a {@code ThreadLocal} so we can use the singleton cache
   * @param readBasesLeft number of bases to left of {@code read()} on read (i.e. to the left of the match)
   */
  @Override
  public void setReadBasesLeft(int readBasesLeft) {
    mReadBasesLeft.set(readBasesLeft);
  }

  /**
   * Be very careful when using this field, it is stored in a {@code ThreadLocal} so we can use the singleton cache
   * @param readBasesRight number of bases to right of {@code read()} on read (i.e. to the right of the match)
   */
  @Override
  public void setReadBasesRight(int readBasesRight) {
    mReadBasesRight.set(readBasesRight);
  }

  @Override
  public boolean isForward() {
    return mForward;
  }

  @Override
  public boolean isReadPaired() {
    return mReadPaired;
  }

  @Override
  public boolean isFirst() {
    return mFirst;
  }

  @Override
  public boolean isMated() {
    return mMated;
  }

  @Override
  public boolean isUnmapped() {
    return mUnmapped;
  }

  @Override
  public boolean integrity() {
    super.integrity();
    Exam.assertTrue(description().valid(mRead));
    Exam.assertTrue(0.0 <= mQ && mQ < 1.0 && !Double.isNaN(mQ));
    Exam.assertTrue(0.0 <= mQt && mQt < 1.0 && !Double.isNaN(mQt));
    Exam.assertTrue(0.0 <= mQm && mQm < 1.0 && !Double.isNaN(mQm));
    Exam.assertTrue(mReadPaired || !mMated);
    return true;
  }
}
