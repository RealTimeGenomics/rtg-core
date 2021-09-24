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
package com.rtg.position.output;

import com.rtg.mode.Frame;


/**
 */
public interface GapScorer {

  /**
   * Provides query id.
   * @param queryId id of current query as from {@link com.rtg.reader.SequencesReader}
   * @param frame frame of query.
   */
  void setQueryId(int queryId, Frame frame);

  /**
   * Provides internal build id
   * @param internalId id of current build sequence encoded for frame
   */
  void setInternalBuildId(int internalId);

  /**
   * Return the score for the gap. Bigger = better
   * Note: in the reverse complement case, the read positions are reported as if
   * the read is reversed, and the template positions as if the template is forward.
   * @param buildStart the start of the region we wish to score on the build (incl)
   * @param buildEnd the end of the region we wish to score on the build (excl)
   * @param queryStart the start of the region we wish to score on the query (incl)
   * @param queryEnd the end of the region we wish to score on the query (excl)
   * @return the score
   */
  double score(final int buildStart, final int buildEnd, final int queryStart, final int queryEnd);

  /**
   * Return the score for the gaps at the ends. Bigger = better
   * Note: in the reverse complement case, the read positions are reported as if
   * the read is reversed, and the template positions as if the template is forward.
   * @param buildStart the start of the region we wish to score on the build (incl)
   * @param buildEnd the end of the region we wish to score on the build (excl)
   * @param queryStart the start of the region we wish to score on the query (incl)
   * @param queryEnd the end of the region we wish to score on the query (excl)
   * @return a maximum score
   */
  double scoreMax(final int buildStart, final int buildEnd, int queryStart, int queryEnd);

  /**
   * Find the maximum delta value that will give a non-zero score.
   * @return the maximum delta (likely to be &gt; 0 but this isn't guaranteed).
   */
  int maxDelta();

  /**
   * Find the minimum delta value that will give a non-zero score.
   * @return the minimum delta (likely to be &gt; 0 but this isn't guaranteed).
   */
  int minDelta();

  /**
   * Perform any end of life processing (including logging statistics)
   */
  void close();
}
