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

/**
 * Builder for {@link GtQualComplexMultiplierModel}.
 * Ignores the VCF files passed to the build method and uses multiplier settings that are given as model parameters.
 *
 */
public class GtQualComplexMultiplierModelBuilder extends AbstractModelBuilder<GtQualComplexMultiplierModel> {

  /** Default GQ score multiplier for simple homozygous */
  public static final double MULTIPLIER_GQ_SIMPLE_HOMOZYGOUS = 2.0;
  /** Default GQ score multiplier for simple heterozygous */
  public static final double MULTIPLIER_GQ_SIMPLE_HETEROZYGOUS = 0.512;
  /** Default GQ score multiplier for complex homozygous */
  public static final double MULTIPLIER_GQ_COMPLEX_HOMOZYGOUS = 1.16;
  /** Default GQ score multiplier for complex heterozygous */
  public static final double MULTIPLIER_GQ_COMPLEX_HETEROZYGOUS = 0.242;

  /** Default QUAL score multiplier for simple calls */
  public static final double MULTIPLIER_QUAL_SIMPLE = 0.2;
  /** Default QUAL score multiplier for complex calls */
  public static final double MULTIPLIER_QUAL_COMPLEX = 0.04;

  /** GQ simple homozygous multiplier parameter property */
  public static final String PARAMETER_MULTIPLIER_SIMPLE_HOMOZYGOUS = "multiplier.gq.simple.homozygous";
  /** GQ simple heterozygous multiplier parameter property */
  public static final String PARAMETER_MULTIPLIER_SIMPLE_HETEROZYGOUS = "multiplier.gq.simple.heterozygous";
  /** GQ complex homozygous multiplier parameter property */
  public static final String PARAMETER_MULTIPLIER_COMPLEX_HOMOZYGOUS = "multiplier.gq.complex.homozygous";
  /** GQ complex heterozygous multiplier parameter property */
  public static final String PARAMETER_MULTIPLIER_COMPLEX_HETEROZYGOUS = "multiplier.gq.complex.heterozygous";

  /** QUAL simple multiplier parameter property */
  public static final String PARAMETER_MULTIPLIER_QUAL_SIMPLE = "multiplier.qual.simple";
  /** QUAL complex multiplier parameter property */
  public static final String PARAMETER_MULTIPLIER_QUAL_COMPLEX = "multiplier.qual.complex";

  /**
   * Create a genotype quality complex call multiplier model.  Ignores attributes and uses <code>GT</code> and <code>XRX</code> fields.
   * Multiplier values are obtained through parameter settings.
   *
   * @param formatAttributes ignored
   * @param infoAttributes ignored
   * @param derivedAttributes ignored
   */
  public GtQualComplexMultiplierModelBuilder(String[] formatAttributes, String[] infoAttributes, String[] derivedAttributes) {
    super(formatAttributes, infoAttributes, derivedAttributes);
    mProperties.setProperty(MODEL_PROPERTY_TYPE, ModelType.GT_COMPLEX.toString());
  }

  @Override
  public void build(VcfDataset... vcfDatasets) {
    // ignores the input vcf files and gets values from parameters

    final GtQualComplexMultiplierModel model = new GtQualComplexMultiplierModel();

    model.setGqSimpleHomozygousMultiplier(getPropertyValue(PARAMETER_MULTIPLIER_SIMPLE_HOMOZYGOUS, MULTIPLIER_GQ_SIMPLE_HOMOZYGOUS));
    model.setGqSimpleHeterozygousMultiplier(getPropertyValue(PARAMETER_MULTIPLIER_SIMPLE_HETEROZYGOUS, MULTIPLIER_GQ_SIMPLE_HETEROZYGOUS));
    model.setGqComplexHomozygousMultiplier(getPropertyValue(PARAMETER_MULTIPLIER_COMPLEX_HOMOZYGOUS, MULTIPLIER_GQ_COMPLEX_HOMOZYGOUS));
    model.setGqComplexHeterozygousMultiplier(getPropertyValue(PARAMETER_MULTIPLIER_COMPLEX_HETEROZYGOUS, MULTIPLIER_GQ_COMPLEX_HETEROZYGOUS));

    model.setQualSimpleMultiplier(getPropertyValue(PARAMETER_MULTIPLIER_QUAL_SIMPLE, MULTIPLIER_QUAL_SIMPLE));
    model.setQualComplexMultiplier(getPropertyValue(PARAMETER_MULTIPLIER_QUAL_COMPLEX, MULTIPLIER_QUAL_COMPLEX));

    mModel = model;
  }

  private double getPropertyValue(String property, double defaultValue) {
    return Double.parseDouble(mParameters.getProperty(property, Double.toString(defaultValue)));
  }

}
