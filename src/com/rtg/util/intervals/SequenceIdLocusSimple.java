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
 * Implementation of a SequenceIdLocus
 */
public class SequenceIdLocusSimple extends Range implements SequenceIdLocus {

  private final int mSequenceId;

  /**
   * @param sequenceId reference sequence identifier
   * @param start position on reference sequence (0 based)
   */
  public SequenceIdLocusSimple(int sequenceId, int start) {
    this(sequenceId, start, start + 1);
  }
  /**
   * @param sequenceId reference sequence identifier
   * @param start position on reference sequence (0 based)
   * @param end position on reference sequence (0 based, exclusive)
   */
  public SequenceIdLocusSimple(int sequenceId, int start, int end) {
    super(start, end);
    mSequenceId = sequenceId;
  }

  @Override
  public int getSequenceId() {
    return mSequenceId;
  }
}
