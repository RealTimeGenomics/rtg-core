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

import java.io.IOException;
import java.util.Iterator;

/**
 * A cache for records fetched from mapping files for one sequence.
 * Records may be fetched more than once except when flush has specified that they will not be requested again.
 *
 * @param <R> type of record to return
 */
public interface ReaderWindow<R extends ReaderRecord<R>> {

  /**
   * Fetch all records whose extent overlaps the specified region. The extent is
   * determined by unrolling the cigar. It is permissible to include records that don't
   * actually overlap if this makes the internal calculations simpler. However, it is expected
   * that these will be kept to a minimum.
   * Callers may receive over-coverage records that need to be dealt with appropriately.
   * @param start first position in region selected (0 based, inclusive).
   * @param end last position in region selected (0 based, exclusive).
   * @return all records selected in any order.
   * @throws IOException if an IO error occurs
   */
  Iterator<R> recordsOverlap(int start, int end) throws IOException;

  /**
   * Explicitly request that records be processed up to the specified point. Ordinarily this happens
   * automatically due to calls to <code>recordsOverlap</code>. But sometimes we want to skip over regions without
   * looking at the records.
   * @param end the end position
   */
  void advanceBuffer(int end);

  /**
   * Specify that all records whose start position is within the specified region
   * will not be requested again by a <code>recordsAtStart</code> or a <code>recordsOverlap</code> call.
   * @param start first position in region selected (0 based, inclusive).
   * @param end last position in region selected (0 based, exclusive).
   * @throws IOException if an IO error occurs
   */
  void flush(int start, int end) throws IOException;

}
