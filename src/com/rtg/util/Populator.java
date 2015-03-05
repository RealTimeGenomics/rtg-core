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
package com.rtg.util;

import htsjdk.samtools.SAMRecord;

/**
 * Interface to populate an object based on the contents of another object.
 * @param <T> target type
 */
public interface Populator<T> {

  /**
   * Construct and populate and object of type <code>T</code> based on the
   * contents of the source object.
   *
   * @param source source object
   * @return populated object
   */
  T populate(SAMRecord source);

  /**
   * A record that can be used to denote an overflow situation.
   * @param position 0-based start position of overflow situation
   * @param length length of overflow region
   * @return overflow indicator
   */
  T overflow(int position, int length);
}

