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

package com.rtg.sam;

import com.rtg.util.intervals.SequenceIdLocus;

/**
 * @param <T> concrete type
 */
public interface ReaderRecord<T extends ReaderRecord<T>> extends SequenceIdLocus, Comparable<T> {

  /**
   * @return which genome the record is from (in range 0 ... ).
   */
  int getGenome();

  /**
   * @return next record in a chain or records or null if there is no more records
   */
  T chain();

  /**
   * @param rec next record to set in the chain
   */
  void setNextInChain(final T rec);

  /**
   * Used to create a consistent ordering on duplicate records.
   * See Bug #1431
   * @param rec the record to compare with
   * @return negative if this record is less than other, zero if equal and positive if greater
   */
  int disambiguateDuplicate(final T rec);
}
