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
package com.rtg.metagenomics.metasnp;

import java.io.IOException;

import com.rtg.mode.DNA;
import com.rtg.util.StringUtils;

/**
 * Represent a variant and allele counts.
 */
public class MetaSnpLine {

  private final String mSequence;
  private final int mPosition;
  final int mReference;
  final String[] mAlleles;
  final double[][] mCounts; // alleles x samples

  private MetaSnpLine(String sequence, int position, int reference, final String[] alleles, final double[][] counts) throws IOException {
    mSequence = sequence;
    mPosition = position;
    mReference = reference;
    mAlleles = alleles;
    mCounts = counts;
  }

  public String getSequence() {
    return mSequence;
  }

  public int getPosition() {
    return mPosition;
  }

  private static final int NUM_FIELDS = 7;
  private static final String[] DEFAULT_ALLELES = {"a", "c", "g", "t"};
  static MetaSnpLine create(final String str, final int lineNumber) throws IOException {
    final String[] fields = StringUtils.split(str, '\t');
    if (fields.length < NUM_FIELDS) {
      throw new IOException("Error on line: " + lineNumber + " Expected at least " + NUM_FIELDS + " columns but found " + fields.length + " ");
    }
    try {
      final double[][] counts = new double[4][];
      for (int i = 0; i < 4; i++) {
        final String count = fields[3 + i];
        final String[] samples = StringUtils.split(count, ',');
        counts[i] = new double[samples.length];
        for (int j = 0; j < samples.length; j++) {
          counts[i][j] = Double.parseDouble(samples[j]);
        }
      }
      return new MetaSnpLine(fields[0], Integer.parseInt(fields[1]) - 1, (byte) DNA.valueOf(fields[2]).ordinal() - 1, DEFAULT_ALLELES, counts);
    } catch (final NumberFormatException e) {
      throw new IOException("Error on line: " + lineNumber + " " + e.getMessage());
    }
  }
}
