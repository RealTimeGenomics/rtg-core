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

/**
 * Information relating to a variant sample
 */
public class VariantSampleInfo {

  /**
   * Enumeration of variant sample format Info possible types
   */
  public enum VariantFormatEnum {
    /** The coverage for a variant */
    COVERAGE,
    /** The coverage correction for a variant */
    COVERAGE_CORRECTION,
    /** The ambiguity ratio for a variant */
    AMBIGUITY_RATIO,
    /** The non identity posterior for a variant */
    NONIDENTITY_POSTERIOR,
    /** The statistics for a variant */
    STATISTICS,
    /** Multiplier that was applied to GQ */
    /** bound on probability that strand bias exists, for allele 1 */
    HOEFFDING_STRAND_BIAS_A1,
    /** bound on probability that strand bias exists, for allele 2 */
    HOEFFDING_STRAND_BIAS_A2,
    /** Hoeffding allele balance for heterozygous calls, bound to probability 0.5 */
    HOEFFDING_ALLELE_BALANCE_HET,
    /** Hoeffding allele balance for homozygous calls, bound to probability 1.0 */
    HOEFFDING_ALLELE_BALANCE_HOM,
    /** bound on probability that site isn't covered evenly */
    HOEFFDING_READ_POSITION,
    /** bound on probability from total unmated/mated ratio, for allele 1 */
    HOEFFDING_UNMATED_BIAS_A1,
    /** bound on probability from total unmated/mated ratio, for allele 2 */
    HOEFFDING_UNMATED_BIAS_A2,
    /** Ratio of placed unmapped reads to mapped reads */
    UNPLACED_RATIO
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
