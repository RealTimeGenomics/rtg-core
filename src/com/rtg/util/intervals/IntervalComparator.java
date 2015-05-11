/*
 * Copyright (c) 2015. Real Time Genomics Limited.
 *
 * Use of this source code is bound by the Real Time Genomics Limited Software Licence Agreement
 * for Academic Non-commercial Research Purposes only.
 *
 * If you did not receive a license accompanying this file, a copy must first be obtained by email
 * from support@realtimegenomics.com.  On downloading, using and/or continuing to use this source
 * code you accept the terms of that license agreement and any amendments to those terms that may
 * be made from time to time by Real Time Genomics Limited.
 */
package com.rtg.util.intervals;

import java.io.Serializable;
import java.util.Comparator;

/**
 * This comparator only uses start and end positions, so should only be used in cases where it is known that all objects
 * being compared reside on the same genomic sequence.
 */
public final class IntervalComparator implements Comparator<Interval>, Serializable {

  @Override
  public int compare(Interval o1, Interval o2) {

    if (o1.getStart() < o2.getStart()) {
      return -1;
    } else if (o1.getStart() > o2.getStart()) {
      return 1;
    }
    if (o1.getEnd() < o2.getEnd()) {
      return -1;
    } else if (o1.getEnd() > o2.getEnd()) {
      return 1;
    }
    return 0;
  }
}
