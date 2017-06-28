/*
 * Copyright (c) 2017. Real Time Genomics Limited.
 *
 * Use of this source code is bound by the Real Time Genomics Limited Software Licence Agreement
 * for Academic Non-commercial Research Purposes only.
 *
 * If you did not receive a license accompanying this file, a copy must first be obtained by email
 * from support@realtimegenomics.com.  On downloading, using and/or continuing to use this source
 * code you accept the terms of that license agreement and any amendments to those terms that may
 * be made from time to time by Real Time Genomics Limited.
 */
package com.rtg.variant.bayes.complex;

import java.util.Collections;
import java.util.List;

import com.rtg.variant.Variant;
import com.rtg.variant.bayes.multisample.VariantAlleleTrigger;

/**
 * A decomposer that only does trimming.
 */
public final class Trimmer extends Splitter {

  /**
   * Construct a new trimmer.
   * @param variantAlleleTrigger method for handling variant alleles
   */
  public Trimmer(VariantAlleleTrigger variantAlleleTrigger) {
    super(null, variantAlleleTrigger);
  }

  @Override
  public List<Variant> decompose(final Variant original) {
    return Collections.singletonList(trim(original, mVariantAlleleTrigger));
  }
}
