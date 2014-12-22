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

import net.sf.samtools.SAMRecord;

/**
 * It's a chain of SAM filters
 */
public class SamFilterChain implements SamFilter {

  private final SamFilter[] mFilters;

  /**
   * Construct a chain of SAM filters
   * @param filters the filters to chain
   */
  public SamFilterChain(SamFilter... filters) {
    mFilters = filters;
  }

  @Override
  public boolean acceptRecord(SAMRecord rec) {
    for (SamFilter filter : mFilters) {
      if (!filter.acceptRecord(rec)) {
        return false;
      }
    }
    return true;
  }
}
