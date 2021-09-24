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
package com.rtg.variant.cnv.region;

import com.rtg.util.intervals.SequenceIdLocusSimple;

/**
 */
abstract class AbstractCnvRegion extends SequenceIdLocusSimple implements Comparable<AbstractCnvRegion>, Region {

  private static int sID;
  private static synchronized int nextID() {
    return ++sID;
  }

  private final int mId = nextID();

  /**
   * @param start position (inclusive)
   * @param end position (inclusive)
   */
  AbstractCnvRegion(final int start, final int end) {
    super(nextID(), start, end);
  }

  @Override
  public int compareTo(final AbstractCnvRegion that) {
    if (that == this) {
      return 0;
    }
    if (this.getStart() < that.getStart()) {
      return -1;
    }
    if (this.getStart() > that.getStart()) {
      return +1;
    }
    //break symmetry in a portable repeatable way

    final int ediff = this.mId - that.mId;
    assert ediff != 0;
    return ediff;
  }


  @Override
  public boolean equals(final Object obj) {
    return super.equals(obj);
  }

  @Override
  public int hashCode() {
    return super.hashCode();
  }



}
