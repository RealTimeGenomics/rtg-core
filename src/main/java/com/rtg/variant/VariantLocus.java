/*
 * Copyright (c) 2018. Real Time Genomics Limited.
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
