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

import java.io.Closeable;
import java.util.Iterator;

import htsjdk.samtools.SAMFileHeader;

/**
 * @param <T> record type
 */
public interface RecordIterator<T> extends Iterator<T>, Closeable {

  /**
   * Get the header for checking.
   * @return the header.
   */
  SAMFileHeader header();

  /**
   * Gets the total number of records that were invalid.
   * @return the sum of all invalid counts
   */
  long getInvalidRecordsCount();

  /**
   * Gets the number of records that were ignored due to filtering criteria
   * @return the count of records ignored due to user-filtering
   */
  long getFilteredRecordsCount();

  /**
   * Gets the number of records that were detected as duplicates and ignored
   * @return the number of duplicate records filtered
   */
  long getDuplicateRecordsCount();

  /**
   * Gets the total number of records that were returned to the caller.
   * @return the count of records returned to the caller.
   */
  long getOutputRecordsCount();

  /**
   * Gets the total number of input records.
   * @return the count of all input records (regardless of validity or filtering status).
   */
  long getTotalRecordsCount();

  /**
   * Get the total number of nucleotides read (ignoring badly formatted records).
   * @return the sum of all nucleotides.
   */
  long getTotalNucleotides();

}
