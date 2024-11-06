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

import static com.rtg.mode.TranslatedFrame.NUCLEOTIDES_PER_CODON;

import com.rtg.index.hash.ngs.HashFunctionFactory;
import com.rtg.index.hash.ngs.general.Skeleton;
import com.rtg.index.hash.ngs.protein.ProteinMask;

/**
 * Mask parameters for protein masks
 */
public class NgsMaskParamsProtein extends NgsMaskParamsGeneral {

  private final boolean mTranslated;

  /**
   * as in super class
   * @param wordSize word size
   * @param substitutions number of substitutions
   * @param indels number of indels
   * @param indelLength length of indels
   * @param translated true for translated (DNA) search, false for untranslated (protein) search
   */
  public NgsMaskParamsProtein(final int wordSize, final int substitutions, final int indels, final int indelLength, boolean translated) {
    super(wordSize, substitutions, indels, indelLength);
    mTranslated = translated;
  }

  @Override
  public HashFunctionFactory maskFactory(int readLength) {
    final Skeleton sk = new Skeleton(convertedReadLength(readLength), mWordSize, mSubstitutions, mIndels, mIndelLength);
    return ProteinMask.factory(sk);
  }

  @Override
  public boolean isValid(int readLength) {
    return super.isValid(convertedReadLength(readLength));
  }

  private int convertedReadLength(int readLength) {
    return (mTranslated ? (readLength / NUCLEOTIDES_PER_CODON) : readLength) - 1;
  }
}
