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

import com.rtg.util.intervals.RegionRestriction;

import htsjdk.samtools.SAMSequenceDictionary;

/**
 * Class to handle SAM region restrictions
 */
public class SamRegionRestriction extends RegionRestriction {

  /**
   * Restriction to named template sequence, accepts restrictions of
   * the form name, name:start-end or name:start+length.
   * @param restriction restriction string
   */
  public SamRegionRestriction(String restriction) {
    super(restriction);
  }

  /**
   * Restriction to named template sequence, accepts restrictions of
   * the form name, name:start-end or name:start+length.
   * @param restriction restriction to be converted
   */
  public SamRegionRestriction(RegionRestriction restriction) {
    super(restriction.getSequenceName(), restriction.getStart(), restriction.getEnd());
  }

  /**
   * Restriction to region in named template
   * @param template sequence name
   * @param start start position of allowed region (0 based inclusive)
   * @param end end position of allowed region (0 based exclusive)
   */
  public SamRegionRestriction(String template, int start, int end) {
    super(template, start, end);
  }

  /**
   * Resolve any missing start position to the actual start position.
   * @return start position for restriction. 0 based inclusive
   */
  public int resolveStart() {
    return (getStart() == MISSING) ? 0 : getStart();
  }

  /**
   * Resolve any missing end position to the actual end position.
   * @param dict end position for restriction resolved to actual position
   * if a wildcard had been used.
   * @return end position for restriction. 0 base exclusive
   */
  public int resolveEnd(SAMSequenceDictionary dict) {
    final int sequenceLength = dict.getSequence(getSequenceName()).getSequenceLength();
    return (getEnd() == MISSING) ? sequenceLength : getEnd();
  }
}
