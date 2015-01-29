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
