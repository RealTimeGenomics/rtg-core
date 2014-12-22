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
package com.rtg.calibrate;

/**
 * Keeps track of statistics about <code>SNP</code> and <code>Indel</code>
 * rates.
 */
public class CalibrationStats {
  private final int[] mValues;
  private long mEqual, mDifferent, mInserted, mDeleted;

  /**
   * @param values the index values for this statistics, can be null for accumulators
   */
  public CalibrationStats(int[] values) {
    mValues = values;
  }

  /**
   * Get the covariate value for the given covariate index.
   * Can only use if non-null value array used to construct.
   * @param covariateIndex the index of the covariate to get the value for
   * @return the value for the given covariate index
   */
  public int getCovariateValue(int covariateIndex) {
    assert mValues != null && covariateIndex < mValues.length;
    return mValues[covariateIndex];
  }

  /**
   * Get the value array that this <code>CalibrationStats</code> represents
   * @return the value array
   */
  int[] getValues() {
    return mValues;
  }

  /**
   * Get the new array index for this this <code>CalibrationStats</code> object using given covariates.
   * Can only use if non-null value array used to construct.
   * @param covariates the covariates to calculate the index for this <code>CalibrationStats</code>
   * @return the index for the given covariates of the values this <code>CalibrationStats</code> represents
   */
  public int getNewIndex(Covariate[] covariates) {
    assert covariates != null && mValues != null && covariates.length == mValues.length;
    int pos = 0;
    for (int i = 0; i < covariates.length; i++) {
      pos = pos * covariates[i].newSize() + mValues[i];
    }
    return pos;
  }

  /**
   * @param length length of the equals region.
   */
  public void seenEquals(final long length) {
    mEqual += length;
  }

  /**
   * @param length length of the SNP/MNP region.
   */
  public void seenMnp(final long length) {
    mDifferent += length;
  }

  /**
   * @param length length of the Delete region.
   */
  public void seenDelete(final long length) {
    mDeleted += length;
  }

  /**
   * @param length length of the Insert region.
   */
  public void seenInsert(final long length) {
    mInserted += length;
  }

  /**
   * @return sum of all stats
   */
  public long getTotalLength() {
    return mEqual + mDifferent + mInserted + mDeleted;
  }

  /**
   * add additional stats into this one
   * @param s stats to add
   */
  public void accumulate(CalibrationStats s) {
    mEqual += s.mEqual;
    mDeleted += s.mDeleted;
    mDifferent += s.mDifferent;
    mInserted += s.mInserted;
  }

  /**
   * @return count deleted
   */
  public long getDeleted() {
    return mDeleted;
  }

  /**
   * @return count different
   */
  public long getDifferent() {
    return mDifferent;
  }

  /**
   * @return count equal
   */
  public long getEqual() {
    return mEqual;
  }

  /**
   * @return count inserted
   */
  public long getInserted() {
    return mInserted;
  }

  /**
   * Convert the values that this <code>CalibrationStats</code> represents back to an index.
   * Can only use if non-null value array used to construct.
   * @param covariates the covariates to use to produce the name with
   * @return the name of the <code>CalibrationStats</code>
   */
  public String getName(Covariate[] covariates) {
    assert covariates != null && mValues != null && covariates.length == mValues.length;
    String sep = "";
    final StringBuilder sb = new StringBuilder();
    for (int i = 0; i < covariates.length; i++) {
      sb.append(sep);
      sb.append(covariates[i].valueString(mValues[i]));
      sep = "\t";
    }
    return sb.toString();
  }

  /**
   * Convert this <code>CalibrationStats</code> to its output string.
   * Can only use if non-null value array used to construct.
   * @param covariates the covariates to use to produce the name with
   * @return the output string for this <code>CalibrationStats</code>
   */
  public String outputString(Covariate[] covariates) {
    return getName(covariates) + "\t" + mEqual + "\t" + mDifferent + "\t" + mInserted + "\t" + mDeleted;
  }
}

