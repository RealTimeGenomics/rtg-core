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

package com.rtg.variant.avr;

import static com.rtg.util.StringUtils.LS;

import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.util.Arrays;
import java.util.Properties;

import com.rtg.util.MathUtils;
import com.rtg.util.Utils;
import com.rtg.variant.format.VcfInfoField;
import com.rtg.vcf.VcfRecord;
import com.rtg.vcf.VcfUtils;
import com.rtg.vcf.header.VcfHeader;

/**
 * Generates a rescaled GQ and QUAL scores based on GT value and whether call is simple or complex.
 * This model ignores the field setting as it updates the GQ and QUAL values.
 */
public final class GtQualComplexMultiplierModel extends AbstractPredictModel {

  private static final String PROP_MULTIPLIER_GQ_SIMPLE_HOMOZYGOUS = "multiplier.gq.simple.homozygous";
  private static final String PROP_MULTIPLIER_GQ_SIMPLE_HETEROZYGOUS = "multiplier.gq.simple.heterozygous";
  private static final String PROP_MULTIPLIER_GQ_COMPLEX_HOMOZYGOUS = "multiplier.gq.complex.homozygous";
  private static final String PROP_MULTIPLIER_GQ_COMPLEX_HETEROZYGOUS = "multiplier.gq.complex.heterozygous";

  private static final String PROP_MULTIPLIER_QUAL_SIMPLE = "multiplier.qual.simple";
  private static final String PROP_MULTIPLIER_QUAL_COMPLEX = "multiplier.qual.complex";

  private static final String PROPERTIES_COMMENT = "QUAL/GQ/complex region multiplier model.";

  // genotype/complexity multipliers
  // order is: simplehomo, simplehet, complexhomo, complexhet
  private final double[] mGqMultipliers = new double[4]; // {1.0, 1.0, 1.0, 1.0};
  private final double[] mQualMultipliers = new double[2]; // {1.0, 1.0};

  private final long[] mSummaryCounts = new long[6];

  /**
   * Default constructor. All multipliers set to 1.0.
   */
  public GtQualComplexMultiplierModel() {
    Arrays.fill(mGqMultipliers, 1.0);
    Arrays.fill(mQualMultipliers, 1.0);
  }

  /**
   * Creates a model getting its settings from the model file.
   *
   * @param is input stream to read from
   * @throws IOException if an error occurs loading from stream
   */
  public GtQualComplexMultiplierModel(InputStream is) throws IOException {
    load(is);
  }

  private void load(InputStream is) throws IOException {
    final Properties props = new Properties();
    props.load(is);
    setGqSimpleHomozygousMultiplier(Double.parseDouble(props.getProperty(PROP_MULTIPLIER_GQ_SIMPLE_HOMOZYGOUS)));
    setGqSimpleHeterozygousMultiplier(Double.parseDouble(props.getProperty(PROP_MULTIPLIER_GQ_SIMPLE_HETEROZYGOUS)));
    setGqComplexHomozygousMultiplier(Double.parseDouble(props.getProperty(PROP_MULTIPLIER_GQ_COMPLEX_HOMOZYGOUS)));
    setGqComplexHeterozygousMultiplier(Double.parseDouble(props.getProperty(PROP_MULTIPLIER_GQ_COMPLEX_HETEROZYGOUS)));
    setQualSimpleMultiplier(Double.parseDouble(props.getProperty(PROP_MULTIPLIER_QUAL_SIMPLE)));
    setQualComplexMultiplier(Double.parseDouble(props.getProperty(PROP_MULTIPLIER_QUAL_COMPLEX)));
  }

  @Override
  public void save(OutputStream os) throws IOException {
    final Properties props = new Properties();
    props.setProperty(PROP_MULTIPLIER_GQ_SIMPLE_HOMOZYGOUS, Double.toString(mGqMultipliers[0]));
    props.setProperty(PROP_MULTIPLIER_GQ_SIMPLE_HETEROZYGOUS, Double.toString(mGqMultipliers[1]));
    props.setProperty(PROP_MULTIPLIER_GQ_COMPLEX_HOMOZYGOUS, Double.toString(mGqMultipliers[2]));
    props.setProperty(PROP_MULTIPLIER_GQ_COMPLEX_HETEROZYGOUS, Double.toString(mGqMultipliers[3]));
    props.setProperty(PROP_MULTIPLIER_QUAL_SIMPLE, Double.toString(mQualMultipliers[0]));
    props.setProperty(PROP_MULTIPLIER_QUAL_COMPLEX, Double.toString(mQualMultipliers[1]));
    props.store(os, PROPERTIES_COMMENT);
  }

  protected void setGqSimpleHomozygousMultiplier(double value) {
    //System.err.println("Setting simple homozygous multiplier to " + value);
    mGqMultipliers[0] = value;
  }

  protected void setGqSimpleHeterozygousMultiplier(double value) {
    //System.err.println("Setting simple heterozygous multiplier to " + value);
    mGqMultipliers[1] = value;
  }

  protected void setGqComplexHomozygousMultiplier(double value) {
    //System.err.println("Setting complex homozygous multiplier to " + value);
    mGqMultipliers[2] = value;
  }

  protected void setGqComplexHeterozygousMultiplier(double value) {
    //System.err.println("Setting complex heterozygous multiplier to " + value);
    mGqMultipliers[3] = value;
  }

  protected void setQualSimpleMultiplier(double value) {
    //System.err.println("Setting QUAL simple multiplier to " + value);
    mQualMultipliers[0] = value;
  }

  protected void setQualComplexMultiplier(double value) {
    //System.err.println("Setting QUAL complex multiplier to " + value);
    mQualMultipliers[1] = value;
  }

  @Override
  public void annotate(VcfRecord record) {
    // is record complex
    final boolean complex = record.hasInfo(VcfInfoField.XRX.name());

    // update QUAL value based on complexity
    final int qualIndex = complex ? 1 : 0;
    final double newQual = Double.parseDouble(record.getQuality()) * mQualMultipliers[qualIndex];
    record.setQuality(Utils.realFormat(newQual, 1));
    incrementCount(qualIndex + mGqMultipliers.length);

    // step over each sample updating GQ value
    for (int s = 0; s < record.getNumberOfSamples(); ++s) {
      // is sample heterozygous
      final boolean heterozygous = VcfUtils.isHeterozygous(record, s);

      // get multiplier based on complex and heterozygous
      final int multiplierIndex = (complex ? 2 : 0) + (heterozygous ? 1 : 0);

      final double gq = VcfUtils.getGenotypeQuality(record, s);

      if (!Double.isNaN(gq)) {
        final String value = Long.toString(MathUtils.round(gq * mGqMultipliers[multiplierIndex]));
        record.getFormat(VcfUtils.FORMAT_GENOTYPE_QUALITY).set(s, value);
      }

      // update summary statistics
      incrementCount(multiplierIndex);
    }
  }

  @Override
  public void annotateSample(VcfRecord record, int sampleNo) {
    final boolean complex = record.hasInfo(VcfInfoField.XRX.name());
    // is sample heterozygous
    final boolean heterozygous = VcfUtils.isHeterozygous(record, sampleNo);

    // get multiplier based on complex and heterozygous
    final int multiplierIndex = (complex ? 2 : 0) + (heterozygous ? 1 : 0);

    final double gq = VcfUtils.getGenotypeQuality(record, sampleNo);

    if (!Double.isNaN(gq)) {
      final String value = Long.toString(MathUtils.round(gq * mGqMultipliers[multiplierIndex]));
      record.getFormat(VcfUtils.FORMAT_GENOTYPE_QUALITY).set(sampleNo, value);
    }

    // update summary statistics
    incrementCount(multiplierIndex);
  }

  private synchronized void incrementCount(int index) {
    mSummaryCounts[index]++;
  }

  @Override
  public void updateHeader(VcfHeader header) {
    // Nothing to add as modifies existing QUAL and GQ fields.
    //header.addFormatField(getField(), MetaType.FLOAT, VcfNumber.ONE, "GT/Complex variant score");
  }

  @Override
  public String toString() {
    assert mGqMultipliers.length == 4;
    assert mQualMultipliers.length == 2;
    return "GQ multipliers:" + LS + PROP_MULTIPLIER_GQ_SIMPLE_HOMOZYGOUS + "\t" + mGqMultipliers[0] + LS + PROP_MULTIPLIER_GQ_SIMPLE_HETEROZYGOUS + "\t" + mGqMultipliers[1] + LS + PROP_MULTIPLIER_GQ_COMPLEX_HOMOZYGOUS + "\t" + mGqMultipliers[2] + LS + PROP_MULTIPLIER_GQ_COMPLEX_HETEROZYGOUS + "\t" + mGqMultipliers[3] + LS + "QUAL multipliers:" + LS + PROP_MULTIPLIER_QUAL_SIMPLE + "\t" + mQualMultipliers[0] + LS + PROP_MULTIPLIER_QUAL_COMPLEX + "\t" + mQualMultipliers[1] + LS;
  }

  @Override
  public String getSummary() {
    assert mSummaryCounts.length == 6;
    final long total = mSummaryCounts[0] + mSummaryCounts[1] +  mSummaryCounts[2] + mSummaryCounts[3];
    final long total2 = mSummaryCounts[4] + mSummaryCounts[5];
    return "GQ multiplier Counts:" + LS + "Simple homozygous\t" + mSummaryCounts[0] + LS + "Simple heterozygous\t" + mSummaryCounts[1] + LS + "Complex homozygous\t" + mSummaryCounts[2] + LS + "Complex heterozygous\t" + mSummaryCounts[3] + LS + "Total\t" + total + LS + "QUAL multiplier Counts:" + LS + "Simple\t" + mSummaryCounts[4] + LS + "Complex\t" + mSummaryCounts[5] + LS + "Total\t" + total2 + LS;
  }

}
