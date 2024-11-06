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

package com.rtg.calibrate;

import com.rtg.sam.BadSuperCigarException;

import htsjdk.samtools.SAMRecord;

/**
 * Interface for variables that might influence the quality of a base.
 *
 */
public interface Covariate {

  /**
   * @return The human readable name of this variable.
   */
  String name();

  /**
   * @return The total number of possible values of this variable.
   */
  int size();

  /**
   * Set the size to be the new size.
   * This should be called after the statistics array has been expanded.
   */
  void resized();

  /**
   * If this is greater than <code>size</code> then the statistics array needs expanding.
   * @return The new total number of possible values of this variable.
   */
  int newSize();

  /**
   * If this is <code>true</code> then the statistics array needs expanding.
   * @return <code>true</code> if the total number of possible values of this variable has increased.
   */
  boolean sizeChanged();

  /**
   * @param sam Sam record for the current read.
   * @param parser This says how far we are through the cigar, the read and the template.
   * @return the value of this variable at the given read and position. (<code>0 .. size()-1</code>)
   * @throws BadSuperCigarException may be thrown.
   */
  int value(SAMRecord sam, CalibratorCigarParser parser) throws BadSuperCigarException;

  /**
   * @param value one of the values that this variable can take (<code>0 .. size()-1</code>)
   * @return the corresponding string for that value.
   */
  String valueString(int value);

  /**
   * @param value one of the string values that this variable can take. may add unknown values to set automatically
   * @return the corresponding index value.
   */
  int parse(String value);

  /**
   * @param value one of the string values that this variable can take.
   * @return the corresponding index value.
   */
  int valueOf(String value);

  /**
   * @return the covariate type
   */
  CovariateEnum getType();
}
