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
