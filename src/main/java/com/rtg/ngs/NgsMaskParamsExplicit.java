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
