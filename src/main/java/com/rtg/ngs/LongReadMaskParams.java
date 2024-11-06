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
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.Integrity;

/**
 * Simply stores the values for long read, is incapable of generating a mask.
 */
public class LongReadMaskParams implements NgsMaskParams, Integrity {

  private final int mWordSize;
  private final int mSubstitutions;
  private final int mIndels;
  private final int mIndelLength;

  /**
   * @param wordSize word size for indexing.
   * @param substitutions guaranteed number of substitutions that will be matched.
   * @param indels guaranteed number of indels that will be matched.
   * @param indelLength maximum length guaranteed to be found for each indel.
   */
  public LongReadMaskParams(final int wordSize, final int substitutions, final int indels, final int indelLength) {
    mWordSize = wordSize;
    mSubstitutions = substitutions;
    mIndels = indels;
    mIndelLength = indelLength;
  }

  @Override
  public int getIndelLength() {
    return mIndelLength;
  }

  @Override
  public int getIndels() {
    return mIndels;
  }

  @Override
  public int getSubstitutions() {
    return mSubstitutions;
  }

  @Override
  public int getWordSize() {
    return mWordSize;
  }

  @Override
  public boolean isValid(int readLength) {
    return true;
  }

  @Override
  public HashFunctionFactory maskFactory(int readLength) {
    throw new UnsupportedOperationException("Not supported, ever");
  }

  @Override
  public void close() {
  }


  @Override
  public String toString() {
    return "Long read Mask: w=" + mWordSize + " s=" + mSubstitutions + " i=" + mIndels + " l=" + mIndelLength;
  }

  @Override
  public boolean globalIntegrity() {
    integrity();
    return true;
  }

  @Override
  public boolean integrity() {
    Exam.assertTrue(0 <= mSubstitutions);
    Exam.assertTrue(0 <= mIndels && mIndels <= mSubstitutions);
    Exam.assertTrue(mIndels == 0 || 1 <= mIndelLength);
    return true;
  }
}
