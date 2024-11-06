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

package com.rtg.index;

import java.util.Arrays;

/**
 * Class to contain an ordered sparse histogram of frequency counts.
 */
public final class SparseFrequencyHistogram {

  private int[] mFrequencies;
  private long[] mCounts;
  private int mUsed = 0;
  private int mLast = Integer.MIN_VALUE;

  /**
   * Constructor for a sparse frequency histogram object.
   * @param initialArraySize the initial size for the arrays.
   */
  private SparseFrequencyHistogram(int initialArraySize) {
    mFrequencies = new int[initialArraySize];
    mCounts = new long[initialArraySize];
  }

  /**
   * Constructor for a sparse frequency histogram object.
   */
  public SparseFrequencyHistogram() {
    this(0);
  }

  /* used for tests only */
  protected int arrayLengths() {
    return mFrequencies.length;
  }

  protected static int calculateNewArrayLength(int currentSize, int minimumNeeded) {
    final int newLength = currentSize / 2 * 3;
    if (newLength < 0) {
      return Integer.MAX_VALUE;
    } else if (newLength < minimumNeeded) {
      return minimumNeeded;
    }
    return newLength;
  }

  private void resize() {
    final int newLength = calculateNewArrayLength(mFrequencies.length, mUsed + 1);
    mFrequencies = Arrays.copyOf(mFrequencies, newLength);
    mCounts = Arrays.copyOf(mCounts, newLength);
  }

  /**
   * Method to add a given count / frequency.
   * For all calls beyond the first the frequency
   * must be greater than or equal to the last
   * frequency provided. Multiple calls with
   * the same frequency will add the count to
   * the existing count.
   * @param frequency the frequency
   * @param count the count to add
   */
  public void add(int frequency, long count) {
    if (mUsed == 0 || mLast != frequency) {
      if (mLast > frequency) {
        throw new IllegalArgumentException("Frequencies must be added in ascending order.");
      }
      if (mUsed >= mFrequencies.length) {
        resize();
      }
      mFrequencies[mUsed] = frequency;
      mLast = frequency;
      ++mUsed;
    }
    mCounts[mUsed - 1] += count;
  }

  /**
   * Get the number of frequency entries contained in the histogram.
   * @return the number of entries contained in the histogram.
   */
  public int length() {
    return mUsed;
  }

  private void checkBounds(int index) {
    if (index < 0 || index >= mUsed) {
      throw new ArrayIndexOutOfBoundsException(index);
    }
  }

  /**
   * Get the frequency for a given index position.
   * @param index the index into the frequency array.
   * @return the frequency from the array index.
   */
  public int getFrequency(int index) {
    checkBounds(index);
    return mFrequencies[index];
  }

  /**
   * Get the count for a given index position.
   * @param index the index into the count array.
   * @return the count from the array index.
   */
  public long getCount(int index) {
    checkBounds(index);
    return mCounts[index];
  }

  /**
   * Method to take two ordered sparse frequency histograms and merge them.
   * @param histogramA the first histogram to merge.
   * @param histogramB the second histogram to merge.
   * @return the new sparse frequency histogram with the merged data.
   */
  public static SparseFrequencyHistogram merge(SparseFrequencyHistogram histogramA, SparseFrequencyHistogram histogramB) {
    final SparseFrequencyHistogram histogram = new SparseFrequencyHistogram(Math.max(histogramA.length(), histogramB.length()));
    int a = 0;
    int b = 0;
    while (a < histogramA.length() && b < histogramB.length()) {
      final int frequencyA = histogramA.getFrequency(a);
      final int frequencyB = histogramB.getFrequency(b);
      if (frequencyA <= frequencyB) {
        histogram.add(frequencyA, histogramA.getCount(a));
        ++a;
      } else {
        histogram.add(frequencyB, histogramB.getCount(b));
        ++b;
      }
    }
    while (a < histogramA.length()) {
      histogram.add(histogramA.getFrequency(a), histogramA.getCount(a));
      ++a;
    }
    while (b < histogramB.length()) {
      histogram.add(histogramB.getFrequency(b), histogramB.getCount(b));
      ++b;
    }
    return histogram;
  }

  /**
   * Method to take an unordered set of individual frequencies and produce
   * a sparse frequency histogram from it.
   * Will also sort the <code>freqDist</code> array passed in.
   * @param freqDist the individual frequency values.
   * @param numUsed the number of elements in the <code>freqDist</code> array being used.
   * @return the sparse frequency histogram for the given individual frequencies.
   */
  public static SparseFrequencyHistogram fromIndividualFrequencies(int[] freqDist, int numUsed) {
    final SparseFrequencyHistogram histogram = new SparseFrequencyHistogram();
    Arrays.sort(freqDist, 0, numUsed);
    for (int i = 0; i < numUsed;) {
      final int freq = freqDist[i];
      final int start = i;
      while (i < numUsed && freq == freqDist[i]) {
        ++i;
      }
      histogram.add(freq, i - start);
    }
    return histogram;
  }
}
