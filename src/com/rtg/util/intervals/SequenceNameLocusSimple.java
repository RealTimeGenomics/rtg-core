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

/**
 * Implementation of a SequenceNameLocus
 */
public class SequenceNameLocusSimple extends Range implements SequenceNameLocus {

  protected final String mSequence;

  /**
   * @param sequence reference sequence name
   * @param start position on reference sequence (0 based)
   * @param end position on reference sequence (0 based, exclusive)
   */
  public SequenceNameLocusSimple(String sequence, int start, int end) {
    super(start, end);
    if (sequence == null) {
      throw new NullPointerException();
    }
    mSequence = sequence;
  }

  @Override
  public String getSequenceName() {
    return mSequence;
  }

  @Override
  public boolean overlaps(SequenceNameLocus other) {
    return overlaps(this, other);
  }

  @Override
  public boolean contains(String sequence, int pos) {
    return contains(this, sequence, pos);
  }

  /**
   * Test whether the supplied regions overlap each other
   * @param current the first range
   * @param other the other range
   * @return true if this region overlaps the other
   */
  public static boolean overlaps(SequenceNameLocus current, SequenceNameLocus other) {
    if (other.getStart() < 0 || other.getEnd() < 0
        || current.getStart() < 0 || current.getEnd() < 0) {
      throw new IllegalArgumentException();
    }
    if (!current.getSequenceName().equals(other.getSequenceName())) {
      return false;
    }
    if (other.getStart() < current.getEnd()
        && other.getEnd() > current.getStart()) {
      return true;
    }
    return false;
  }


  /**
   * Test whether the supplied position is within the current regions
   * @param current the current range
   * @param sequence the name of the sequence containing the other position
   * @param pos the position within the sequence
   * @return true if this region overlaps the other
   */
  public static boolean contains(SequenceNameLocus current, String sequence, int pos) {
    if (current.getStart() < 0 || current.getEnd() < 0) {
      throw new IllegalArgumentException();
    }
    return current.getSequenceName().equals(sequence) && pos >= current.getStart() && pos < current.getEnd();
  }

  @Override
  public String toString() {
    return getSequenceName() + ":" + (getStart() + 1) + "-" + getEnd(); // This representation (end is 1 based inclusive) is for consistency with RegionRestriction
  }

}
