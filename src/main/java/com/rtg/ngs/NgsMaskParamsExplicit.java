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

import com.rtg.index.hash.ngs.FactoryUtil;
import com.rtg.index.hash.ngs.HashFunctionFactory;
import com.rtg.util.ObjectParams;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.Integrity;

/**
 */
public final class NgsMaskParamsExplicit extends ObjectParams implements NgsMaskParams, Integrity {

  private final String mMask;

  /**
   * @param mask name of mask class.
   */
  public NgsMaskParamsExplicit(final String mask) {
    mMask = mask;
    mObjects = new Object[] {mMask};
  }

  @Override
  public HashFunctionFactory maskFactory(int readLength) {
    return FactoryUtil.hashFunction(mMask);
  }

  @Override
  public String toString() {
    return "Mask:" + mMask;
  }

  @Override
  public boolean integrity() {
    Exam.assertTrue(mMask != null && !mMask.trim().isEmpty());
    return true;
  }

  @Override
  public boolean globalIntegrity() {
    integrity();
    return true;
  }

  @Override
  public int getIndelLength() {
    throw new UnsupportedOperationException("Not supported ever.");
  }

  @Override
  public int getIndels() {
    throw new UnsupportedOperationException("Not supported ever.");
  }

  @Override
  public int getSubstitutions() {
    throw new UnsupportedOperationException("Not supported ever.");
  }

  @Override
  public int getWordSize() {
    throw new UnsupportedOperationException("Not supported ever.");
  }

  @Override
  public boolean isValid(int readLength) {
    return true;
  }

}
