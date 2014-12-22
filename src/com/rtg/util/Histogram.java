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

package com.rtg.util;

import java.util.Arrays;

/**
 * Class for holding an expandable histogram.
 */
public class Histogram {

  private long[] mHistogram;

  /**
   * Constructor
   */
  public Histogram() {
    mHistogram = new long[0];
  }

  /**
   * Increment the value at position by 1.
   * @param position the zero based position to increment
   */
  public void increment(int position) {
    increment(position, 1);
  }

  /**
   * Increment the value at position by value.
   * @param position the zero based position to increment
   * @param value the value to increment by
   */
  public void increment(int position, long value) {
    assert value >= 0;
    if (position >= mHistogram.length) {
      mHistogram = Arrays.copyOf(mHistogram, position + 1);
    }
    mHistogram[position] += value;
  }

  /**
   * Get the length of the histogram.
   * @return the length of the histogram
   */
  public int getLength() {
    return mHistogram.length;
  }

  /**
   * Get the value at the given position.
   * @param position the zero based position
   * @return the value at the given position
   */
  public long getValue(int position) {
    return mHistogram[position];
  }

  /**
   * Get the value at the given position, returning
   * 0 for values greater than the length of the histogram.
   * @param position the zero based position
   * @return the value at the given position
   */
  public long getValueUnbounded(int position) {
    return position >= mHistogram.length ? 0 : position < 0 ? 0 : mHistogram[position];
  }

  @Override
  public String toString() {
    final StringBuilder sb = new StringBuilder();
    for (long i : mHistogram) {
      sb.append(i).append("\t");
    }
    return sb.toString().trim();
  }

  /**
   * Parse a zero based histogram from tab separated string.
   * @param histStr string to parse and add
   */
  public void addHistogram(String histStr) {
    if (histStr.length() > 0) {
      final String[] values = histStr.split("\t");
      for (int i = values.length - 1; i >= 0; i--) {
        final long val = Long.parseLong(values[i]);
        if (val > 0) {
          increment(i, val);
        }
      }
    }
  }

  /**
   * Merges the contents of another histogram into this one.
   * @param other the other histogram
   */
  public void addHistogram(Histogram other) {
    for (int i = 0; i < other.getLength(); i++) {
      increment(i, other.getValue(i));
    }
  }

  /**
   * Create a distribution from current histogram
   * @return the distribution
   */
  public double[] toDistribution() {
    final double[] ret = new double[getLength()];
    long tot = 0;
    for (int i = 0; i < getLength(); i++) {
      tot += getValue(i);
    }
    for (int i = 0; i < getLength(); i++) {
      ret[i] = (double) getValue(i) / tot;
    }
    return ret;
  }
}
