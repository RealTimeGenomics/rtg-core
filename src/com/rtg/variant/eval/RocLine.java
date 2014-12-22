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

package com.rtg.variant.eval;

import com.rtg.util.Utils;

/**
 * Represents a line of ROC output (i.e. a point on a ROC curve)
 */
class RocLine implements Comparable<RocLine> {
  final String mSequence;
  final int mPos;
  final double mPrimarySortValue;
  final double mWeight;
  final boolean mCorrect;

  RocLine(String sequence, int pos, double primarySortValue, double weight, boolean correct) {
    super();
    mSequence = sequence;
    mPos = pos;
    mPrimarySortValue = primarySortValue;
    mWeight = weight;
    mCorrect = correct;
  }

  @Override
  public int compareTo(RocLine other) {
    final int order = Double.compare(other.mPrimarySortValue, this.mPrimarySortValue);
    if (order != 0) {
      return order;
    }
    if (!this.mCorrect && other.mCorrect) {
      return -1;
    } else if (this.mCorrect && !other.mCorrect) {
      return 1;
    }
    final int sequence = this.mSequence.compareTo(other.mSequence);
    if (sequence != 0) {
      return sequence;
    }
    return this.mPos - other.mPos;
  }

  @Override
  public boolean equals(Object other) {
    if (other == null) {
      return false;
    }
    if (!other.getClass().equals(getClass())) {
      return false;
    }
    final RocLine o = (RocLine) other;
    return compareTo(o) == 0;
  }

  @Override
  public int hashCode() {
    return Utils.hash(new Object[] {mSequence, mPos, mPrimarySortValue, mCorrect});
  }

}
