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
