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
package com.rtg.util.intervals;

import java.io.Serializable;
import java.util.Comparator;

/**
 * Comparator for ordering SequenceNameLocus objects.
 */
public class SequenceNameLocusComparator implements Comparator<SequenceNameLocus>, Serializable {

  @Override
  public int compare(final SequenceNameLocus o1, final SequenceNameLocus o2) {
    final int orderSequence = o1.getSequenceName().compareTo(o2.getSequenceName());
    if (orderSequence != 0) {
      return orderSequence;
    }
    if (o1.getStart() < o2.getStart()) {
      return -1;
    } else if (o1.getStart() > o2.getStart()) {
      return 1;
    }
    return o1.getEnd() - o2.getEnd();
  }
}
