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

import com.rtg.util.IntegerOrPercentage;

import htsjdk.samtools.SAMRecord;

/**
 * Filtering based on a filter parameters.
 */
public class DefaultSamFilter implements SamFilter {

  private final SamFilterParams mFilterParams;

  /**
   * Filter for SAM records.
   *
   * @param params filter parameters
   */
  public DefaultSamFilter(final SamFilterParams params) {
    mFilterParams = params;
  }

  /**
   * Filter a SAM record based on specified filtering parameters. Does not consider position based filters.
   *
   * @param params parameters
   * @param rec record
   * @return true if record should be accepted for processing
   */
  public static boolean acceptRecord(final SamFilterParams params, final SAMRecord rec) {
    final int flags = rec.getFlags();
    if ((flags & params.requireUnsetFlags()) != 0) {
      return false;
    }
    if ((flags & params.requireSetFlags()) != params.requireSetFlags()) {
      return false;
    }
    if (rec.getAlignmentStart() == 0 && params.excludeUnplaced()) {
      return false;
    }

    final boolean mated = (flags & SamBamConstants.SAM_READ_IS_MAPPED_IN_PROPER_PAIR) != 0;
    final boolean unmapped = rec.getReadUnmappedFlag();
    if (params.excludeUnmated() && !mated && !unmapped) {
      return false;
    }

    final int minMapQ = params.minMapQ();
    if (minMapQ >= 0 && rec.getMappingQuality() < minMapQ) {
      return false;
    }

    final int maxNH = params.maxAlignmentCount();
    if (maxNH >= 0) {
      final Integer nh = SamUtils.getNHOrIH(rec);
      if (nh != null && nh > maxNH) {
        return false;
      }
    }

    final IntegerOrPercentage maxAS = mated ? params.maxMatedAlignmentScore() : params.maxUnmatedAlignmentScore();
    if (maxAS != null) {
      final Integer as = rec.getIntegerAttribute(SamUtils.ATTRIBUTE_ALIGNMENT_SCORE);
      if (as != null && as > maxAS.getValue(rec.getReadLength())) {
        return false;
      }
    }
    return true;
  }


  @Override
  public boolean acceptRecord(final SAMRecord rec) {
    return acceptRecord(mFilterParams, rec);
  }
}
