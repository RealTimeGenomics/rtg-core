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
package com.rtg.tabix;

import java.io.IOException;

import com.rtg.util.intervals.ReferenceRanges;
import com.rtg.util.intervals.SequenceNameLocus;

/**
 * Common lookup interface for tabix and bam index.
 */
public interface LocusIndex {

  /**
   * Retrieve the file pointers associated with region specified
   * @param region the region to extract offsets for
   * @return file pointers for records, or null if no records are indexed for given range
   * @throws IOException if an IO exception occurs
   */
  VirtualOffsets getFilePointers(SequenceNameLocus region) throws IOException;

  /**
   * Retrieve the file pointers associated with regions specified
   * @param ranges contains the query regions
   * @return file pointers for records, or null if no records are indexed for given range
   * @throws IOException if an IO exception occurs
   */
  VirtualOffsets getFilePointers(ReferenceRanges<String> ranges) throws IOException;


  /**
   * @return the names of the sequences this index file refers to
   */
  String[] sequenceNames();
}
