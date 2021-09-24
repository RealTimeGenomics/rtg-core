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
package com.rtg.variant.avr;

/**
 * Builder for {@link NullModel}.
 * Ignores the VCF files passed to the build method.
 *
 */
public class NullModelBuilder extends AbstractModelBuilder<NullModel> {

  /**
   * Create a null model.  Ignores everything.
   *
   * @param formatAttributes ignored
   * @param infoAttributes ignored
   * @param derivedAttributes ignored
   */
  public NullModelBuilder(String[] formatAttributes, String[] infoAttributes, String[] derivedAttributes) {
    super(formatAttributes, infoAttributes, derivedAttributes);
    mProperties.setProperty(MODEL_PROPERTY_TYPE, ModelType.NULL.toString());
  }

  @Override
  public void build(VcfDataset... vcfDatasets) {
    // ignores the input vcf files and gets values from parameters
    mModel = new NullModel();
  }
}
