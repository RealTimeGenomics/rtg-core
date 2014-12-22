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
 * Comparator for ordering SequenceIdLocus objects.
 */
public class SequenceIdLocusComparator implements Comparator<SequenceIdLocus>, Serializable {

  @Override
  public int compare(SequenceIdLocus o1, SequenceIdLocus o2) {
    int cmp = o1.getSequenceId() - o2.getSequenceId();
    if (cmp != 0) {
      return cmp;
    }

    cmp = o1.getStart() - o2.getStart();
    if (cmp != 0) {
      return cmp;
    }

    return o1.getEnd() - o2.getEnd();
  }
}

