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

import java.util.List;

/**
 * Container for an entire index to be written
 */
public class SequenceIndexContainer {
  List<SequenceIndex> mIndexes;
  long mNumUnmappedNoCoordinates;

  /**
   * Constructor
   * @param indexes list of indices
   * @param numUnmapped number of unmapped records which had no coordinates
   */
  public SequenceIndexContainer(List<SequenceIndex> indexes, long numUnmapped) {
    mIndexes = indexes;
    mNumUnmappedNoCoordinates = numUnmapped;
  }

  /**
   * @return get indices
   */
  public List<SequenceIndex> getIndexes() {
    return mIndexes;
  }
  /**
   * @return number of unmapped records which had no coordinates
   */
  public long numUnmappedNoCoordinates() {
    return mNumUnmappedNoCoordinates;
  }
}
