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
import com.rtg.index.hash.ngs.general.Skeleton;
import com.rtg.index.hash.ngs.protein.ProteinMask;

/**
 * Mask parameters for protein masks
 */
public class NgsMaskParamsProtein extends NgsMaskParamsGeneral {

  /**
   * as in super class
   * @param wordSize word size
   * @param substitutions number of substitutions
   * @param indels number of indels
   * @param indelLength length of indels
   */
  public NgsMaskParamsProtein(final int wordSize, final int substitutions, final int indels, final int indelLength) {
    super(wordSize, substitutions, indels, indelLength, false);
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


  private static int convertedReadLength(int readLength) {
    return readLength / 3 - 1;
  }
}
