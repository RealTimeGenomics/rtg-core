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

/**
 * Information relating to a variant sample.
 */
public class VariantSampleInfo {

  /**
   * Enumeration of variant sample format Info possible types.
   */
  public enum VariantFormatEnum {
    /** The coverage for a variant. */
    COVERAGE,
    /** The coverage correction for a variant. */
    COVERAGE_CORRECTION,
    /** The ambiguity ratio for a variant. */
    AMBIGUITY_RATIO,
    /** The non-identity posterior for a variant. */
    NONIDENTITY_POSTERIOR,
    /** The statistics for a variant. */
    STATISTICS,
    /** Bound on probability that strand bias exists, for allele 1. */
    HOEFFDING_STRAND_BIAS_A1,
    /** Bound on probability that strand bias exists, for allele 2. */
    HOEFFDING_STRAND_BIAS_A2,
    /** Hoeffding allele balance for heterozygous calls, bound to probability 0.5. */
    HOEFFDING_ALLELE_BALANCE_HET,
    /** Hoeffding allele balance for homozygous calls, bound to probability 1.0. */
    HOEFFDING_ALLELE_BALANCE_HOM,
    /** Bound on probability that site isn't covered evenly. */
    HOEFFDING_READ_POSITION,
    /** Bound on probability from total unmated/mated ratio, for allele 1. */
    HOEFFDING_UNMATED_BIAS_A1,
    /** Bound on probability from total unmated/mated ratio, for allele 2. */
    HOEFFDING_UNMATED_BIAS_A2,
    /** Ratio of placed unmapped reads to mapped reads. */
    UNPLACED_RATIO,
    /** Somatic score. */
    SOMATIC_SCORE
  }

  private Double mDoubleVal;
  private String mStringVal;
  private Integer mIntVal;
  private Boolean mBoolVal;

  /**
   * Create a call info with a double
   * @param value the double value
   */
  public VariantSampleInfo(Double value) {
    mDoubleVal = value;
  }

  /**
   * Create a call info with a String
   * @param value the string value
   */
  public VariantSampleInfo(String value) {
    mStringVal = value;
  }

  /**
   * Create a call info with an Integer
   * @param value the integer value
   */
  public VariantSampleInfo(Integer value) {
    mIntVal = value;
  }

  /**
   * Create a call info with a boolean
   * @param value the boolean value
   */
  public VariantSampleInfo(Boolean value) {
    mBoolVal = value;
  }

  String stringValue() {
    return mStringVal;
  }

  Double doubleValue() {
    return mDoubleVal;
  }

  Integer intValue() {
    return mIntVal;
  }

  /**
   * @return the boolean value of this info. Note returns false even if value hasn't been explicitly set!
   */
  Boolean boolValue() {
    return mBoolVal != null && mBoolVal;
  }
}
