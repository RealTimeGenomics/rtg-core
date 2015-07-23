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

package com.rtg.simulation.genome;

import java.io.IOException;

import com.rtg.reader.SequencesReader;

/**
 * Encapsulate the random selection of sequences
 */
public class SequenceDistribution {
  private final double[] mCumulative;

  /**
   * create a sequence distribution based on the ratios provided in <code>nonCumulative</code>
   * @param nonCumulative the ratios between sequences
   */
  public SequenceDistribution(double[] nonCumulative) {
    mCumulative = new double[nonCumulative.length];
    double sum = 0;
    for (int i = 0; i < nonCumulative.length; i++) {
      sum += nonCumulative[i];
      mCumulative[i] = sum;
    }
  }

  /**
   * Pick a sequence based on where in the distribution the provided double falls
   * @param rand double between 0 and 1
   * @return the sequence id
   */
  public int selectSequence(double rand) {
    int seqId = 0;
    while (mCumulative[seqId] < rand) {
      seqId++;
    }
    return seqId;
  }

  /**
   * Create a distribution using either the provided distribution or a default one based on sequence lengths
   * @param sr the sequences this distribution should cover
   * @param prob the desired distribution or null if default should be used
   * @return an instantiated SequenceDistribution
   * @throws IOException if the reader complains
   */
  public static SequenceDistribution createDistribution(SequencesReader sr, double[] prob) throws IOException {
    if (prob == null) {
      return defaultDistribution(sr);
    } else {
      return new SequenceDistribution(prob);
    }
  }

  /**
   * Create a distribution where the odds of selecting a sequence are proportional to it's length
   * @param sr the sequences this distribution should cover
   * @return an instantiated SequenceDistribution
   * @throws IOException if the reader complains
   */
  public static SequenceDistribution defaultDistribution(SequencesReader sr) throws IOException {
    final int[] lengths = sr.sequenceLengths(0, sr.numberSequences());
    long total = 0;
    for (int length : lengths) {
      total += length;
    }
    final double[] nonCumulative = new double[lengths.length];
    for (int i = 0; i < lengths.length; i++) {
      nonCumulative[i] = (double) lengths[i] / (double) total;
    }
    return new SequenceDistribution(nonCumulative);
  }
}
