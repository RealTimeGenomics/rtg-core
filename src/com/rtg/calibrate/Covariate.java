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
