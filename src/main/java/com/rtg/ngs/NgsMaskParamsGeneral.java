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
package com.rtg.ngs;

import com.rtg.index.hash.ngs.HashFunctionFactory;
import com.rtg.index.hash.ngs.general.Mask;
import com.rtg.index.hash.ngs.general.Skeleton;
import com.rtg.util.ObjectParams;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.Integrity;

/**
 */
public class NgsMaskParamsGeneral extends ObjectParams implements NgsMaskParams, Integrity {
  protected final int mWordSize;

  protected final int mSubstitutions;

  protected final int mIndels;

  protected final int mIndelLength;


  /**
   * @param wordSize word size for indexing.
   * @param substitutions guaranteed number of substitutions that will be matched.
   * @param indels guaranteed number of indels that will be matched.
   * @param indelLength maximum length guaranteed to be found for each indel.
   */
  public NgsMaskParamsGeneral(int wordSize, int substitutions, int indels, int indelLength) {
    mWordSize = wordSize;
    mSubstitutions = substitutions;
    mIndels = indels;
    mIndelLength = indelLength;
    mObjects = new Object[] {mWordSize, mSubstitutions, mIndels, mIndelLength};

  }

  /**
   * Get
   * @return word size
   */
  @Override
  public int getWordSize() {
    return mWordSize;
  }

  /**
   * Get
   * @return number of indels
   */
  @Override
  public int getIndels() {
    return mIndels;
  }

  /**
   * Get
   * @return length of indels
   */
  @Override
  public int getIndelLength() {
    return mIndelLength;
  }

  /**
   * Get
   * @return number of substitutions
   */
  @Override
  public int getSubstitutions() {
    return mSubstitutions;
  }

  @Override
  public HashFunctionFactory maskFactory(int readLength) {
    final Skeleton sk = new Skeleton(readLength, mWordSize, mSubstitutions, mIndels, mIndelLength);
    return Mask.factory(sk);
  }

  @Override
  public boolean isValid(int readLength) {
    final Skeleton sk = new Skeleton(readLength, mWordSize, mSubstitutions, mIndels, mIndelLength);
    return sk.valid();
  }

  @Override
  public String toString() {
    return "General Mask: w=" + mWordSize + " s=" + mSubstitutions + " i=" + mIndels + " l=" + mIndelLength;
  }

  @Override
  public boolean integrity() {
    Exam.assertTrue(0 <= mIndels && mIndels <= mSubstitutions);
    Exam.assertTrue(mIndels == 0 || (1 <= mIndelLength));
    return true;
  }

  @Override
  public boolean globalIntegrity() {
    integrity();
    return true;
  }

}
