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
import com.rtg.util.Params;


/**
 */
public interface NgsMaskParams extends Params {

    /**
     * Get a factory for the <code>HashFunction</code> as specified by the command line parameters.
     * @return a factory for the <code>HashFunction</code>.
     * @param readLength the length of reads that masks will be used on
     */
    HashFunctionFactory maskFactory(int readLength);


  /**
   * Get
   * @return word size
   */
  int getWordSize();

  /**
   * Get
   * @return number of indels
   */
  int getIndels();

  /**
   * Get
   * @return length of indels
   */
  int getIndelLength();

  /**
   * Get
   * @return number of substitutions
   */
  int getSubstitutions();

  /**
   * @return true if this mask is valid
   * @param readLength the length of reads that masks will be used on
   */
  boolean isValid(int readLength);
}
