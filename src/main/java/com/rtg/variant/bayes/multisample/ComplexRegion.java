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
package com.rtg.variant.bayes.multisample;

import com.rtg.util.integrity.Exam;
import com.rtg.util.intervals.SequenceNameLocusSimple;

/**
 * Identifies a complex region
 */
public class ComplexRegion extends SequenceNameLocusSimple {

  /**
   * Describes the type of the region
   */
  public enum RegionType {
    /** lone interesting SNP call */
    INTERESTING("IN", "interesting"),
    /** multiple SNP calls in close proximity or a indel */
    COMPLEX("CX", "complex-called"),
    /** very long complex region */
    HYPER("HX", "hyper-complex"),
    /** very high coverage */
    OVERFLOW("EC", "extreme-coverage"),
    /** NO call due to over coverage */
    OVERCOVERAGE("OC", "complex-over-coverage"),
    /** No call due to no hypotheses */
    NO_HYPOTHESES("NH", "complex-no-hypotheses"),
    /** Complex calling resulted in no variants */
    COMPLEX_NO_VARIANT("CXF", "complex-no-variant"),
    /** No call due to too many hyp */
    TOO_MANY_HYPOTHESES("TMH", "complex-too-many-hypotheses");

    private final String mAbbreviation;
    private final byte[] mDescription;

    RegionType(String abbreviation, String name) {
      mAbbreviation = abbreviation;
      mDescription = name.getBytes();
    }

    String abbreviation() {
      return mAbbreviation;
    }

    /**
     * @return description for the region, used in bed output
     */
    public byte[] description() {
      return  mDescription;
    }
  }

  private RegionType mType;
  private final int mIndelStart;
  private final int mIndelEnd;

  /**
   * create a regular complex region (i.e. not dangling)
   * @param refName name of the sequence
   * @param start start position (0 based, inclusive)
   * @param end end position (0 based, exclusive)
   * @param type the type of the region
   */
  public ComplexRegion(String refName, int start, int end, RegionType type) {
    this(refName, start, end, type, 0, 0);
  }

  /**
   * create a regular complex region (i.e. not dangling)
   * @param refName name of the sequence
   * @param start start position (0 based, inclusive)
   * @param end end position (0 based, exclusive)
   * @param type the type of the region
   * @param indelStart the length of indel associated with the region start
   * @param indelEnd the length of indel associated with the region end
   */
  public ComplexRegion(String refName, int start, int end, RegionType type, int indelStart, int indelEnd) {
    super(refName, start, end);
    mType = type;
    mIndelStart = indelStart;
    mIndelEnd = indelEnd;
  }

  /**
   * @return the type of the region
   */
  public RegionType type() {
    return mType;
  }

  /** @return length of indel at start of region */
  int getIndelStart() {
    return mIndelStart;
  }

  /** @return length of indel at end of region */
  int getIndelEnd() {
    return mIndelEnd;
  }

  @Override
  public String toString() {
    return getSequenceName() + "[" + getStart() + ".." + getEnd() + ")" + mType.abbreviation();
  }

  boolean integrity() {
    Exam.assertTrue(getStart() >= 0);
    Exam.assertTrue(getEnd() >= getStart());
    Exam.assertNotNull(mType);
    if (mType == RegionType.INTERESTING) {
      Exam.assertEquals(getEnd(), getStart() + 1);
    }
    if (mType == RegionType.HYPER) {
      Exam.assertTrue(getEnd() - getStart() > 1);
    }
    return true;
  }

  /**
   * @param type type of the region
   */
  public void setType(RegionType type) {
    mType = type;
  }

}
