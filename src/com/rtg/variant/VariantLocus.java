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

package com.rtg.variant;

import com.rtg.util.intervals.SequenceNameLocusSimple;

/**
 * Stores position information associated with a call.
 */
public class VariantLocus extends SequenceNameLocusSimple {

  private final char mPreviousRefNt;
  private final String mRefNts;

  /**
   * Create a bogus VariantLocus that is not ever intended to actually be output.
   * These are currently used for overlap records, and synthetic variants created in complexities.
   *
   * @param sequence The name of the reference sequence
   * @param start the start position on the reference, 0 based
   * @param end the end position on the reference, 0 based exclusive
   */
  public VariantLocus(String sequence, int start, int end) {
    super(sequence, start, end);
    mRefNts = null;
    mPreviousRefNt = '\0';
  }

  /**
   * @param sequence The name of the reference sequence
   * @param start the start position on the reference, 0 based
   * @param end the end position on the reference, 0 based exclusive
   * @param refNts the reference nucleotides used in this call
   * @param prevRefNt the previous reference nucleotide
   */
  public VariantLocus(String sequence, int start, int end, String refNts, char prevRefNt) {
    super(sequence, start, end);
    assert refNts.length() == (end - start) : "Supplied reference '" + refNts + "' does not match length of interval " + start + "-" + end;
    mRefNts = refNts;
    mPreviousRefNt = prevRefNt;
  }

  public String getRefNts() {
    return mRefNts;
  }

  public char getPreviousRefNt() {
    return mPreviousRefNt;
  }

  @Override
  public String toString() {
    return "start=" + getStart() + " end=" + getEnd()/* + " name=" + getReferenceName()*/;
  }

}
