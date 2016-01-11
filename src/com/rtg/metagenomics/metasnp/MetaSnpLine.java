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
import java.util.List;

import com.rtg.mode.DNA;
import com.rtg.util.StringUtils;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.variant.format.VcfFormatField;
import com.rtg.vcf.VcfRecord;

/**
 * Represent a variant and allele counts.
 */
public final class MetaSnpLine {

  private final String mSequence;
  private final int mPosition;
  private final int mReference;
  final String[] mAlleles;
  final double[][] mCounts; // alleles x samples

  private MetaSnpLine(String sequence, int position, int reference, final String[] alleles, final double[][] counts) {
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

  public int getReferenceIndex() {
    return mReference;
  }

  public String getReferenceAllele() {
    return mAlleles[mReference];
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
      for (int i = 0; i < counts.length; i++) {
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

  private static final String AD = VcfFormatField.AD.name();

  static MetaSnpLine create(final VcfRecord rec) {
    final List<String> alts = rec.getAltCalls();
    final String[] alleles = new String[alts.size() + 1];
    alleles[0] = rec.getRefCall();
    for (int k = 1; k < alleles.length; k++) {
      alleles[k] = alts.get(k - 1);
    }
    final double[][] counts = new double[alleles.length][rec.getNumberOfSamples()];
    for (int j = 0; j < rec.getNumberOfSamples(); j++) {
      final String ade = rec.getSampleString(j, AD); // todo ADE
      if (ade == null) {
        throw new UnsupportedOperationException(); // todo AD fallback?
      }
      final String[] perAllele = StringUtils.split(ade, ',');
      if (perAllele.length != alleles.length) {
        throw new NoTalkbackSlimException("Expected " + alleles.length + " entries in AD(E) field of " + rec);
      }
      for (int k = 0; k < perAllele.length; k++) {
        try {
          counts[k][j] = Double.parseDouble(perAllele[k]);
        } catch (final NumberFormatException e) {
          throw new NoTalkbackSlimException("Error parsing " + perAllele[k] + " as a number in VCF record: " + rec);
        }
      }
    }
    return new MetaSnpLine(rec.getSequenceName(), rec.getStart(), 0, alleles, counts);
  }
}
