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
package com.rtg.metagenomics.metasnp;

import static com.rtg.vcf.VcfUtils.FORMAT_ADE;
import static com.rtg.vcf.VcfUtils.FORMAT_ALLELIC_DEPTH;
import static com.rtg.vcf.VcfUtils.MISSING_FIELD;

import java.io.IOException;
import java.util.List;

import com.rtg.mode.DNA;
import com.rtg.util.StringUtils;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
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
  private static final String[] DEFAULT_ALLELES = {"A", "C", "G", "T"};

  static MetaSnpLine create(final String str, final int lineNumber) throws IOException {
    final String[] fields = StringUtils.split(str, '\t');
    if (fields.length < NUM_FIELDS) {
      throw new IOException("Error on line: " + lineNumber + " Expected at least " + NUM_FIELDS + " columns but found " + fields.length + " ");
    }
    try {
      final double[][] counts = new double[4][];
      for (int i = 0; i < counts.length; ++i) {
        final String count = fields[3 + i];
        final String[] samples = StringUtils.split(count, ',');
        counts[i] = new double[samples.length];
        for (int j = 0; j < samples.length; ++j) {
          counts[i][j] = Double.parseDouble(samples[j]);
        }
      }
      return new MetaSnpLine(fields[0], Integer.parseInt(fields[1]) - 1, (byte) DNA.valueOf(fields[2]).ordinal() - 1, DEFAULT_ALLELES, counts);
    } catch (final NumberFormatException e) {
      throw new IOException("Error on line: " + lineNumber + " " + e.getMessage(), e);
    }
  }

  static MetaSnpLine create(final VcfRecord rec) {
    final List<String> alts = rec.getAltCalls();
    final String[] alleles = new String[alts.size() + 1];
    alleles[0] = rec.getRefCall();
    for (int k = 1; k < alleles.length; ++k) {
      alleles[k] = alts.get(k - 1);
    }
    final double[][] counts = new double[alleles.length][rec.getNumberOfSamples()];
    for (int j = 0; j < rec.getNumberOfSamples(); ++j) {
      String ade = rec.getSampleString(j, FORMAT_ADE);
      if (ade == null) {
        ade = rec.getSampleString(j, FORMAT_ALLELIC_DEPTH);
        if (ade == null) {
          throw new NoTalkbackSlimException("A VCF file with ADE or AD format fields is required.");
        }
      }
      if (MISSING_FIELD.equals(ade)) {
        return null; // This record has no AD or ADE, skip it
      }
      final String[] perAllele = StringUtils.split(ade, ',');
      if (perAllele.length != alleles.length) {
        throw new NoTalkbackSlimException("Expected " + alleles.length + " entries in AD(E) field of " + rec);
      }
      for (int k = 0; k < perAllele.length; ++k) {
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
