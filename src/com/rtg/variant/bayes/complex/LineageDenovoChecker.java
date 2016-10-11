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
package com.rtg.variant.bayes.complex;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import com.rtg.relation.LineageLookup;
import com.rtg.util.StringUtils;
import com.rtg.variant.Variant;
import com.rtg.variant.VariantSample;
import com.rtg.variant.util.VariantUtils;

/**
 * De novo checker that supports cell lineage de novo checking
 */
public class LineageDenovoChecker implements DenovoChecker {

  final LineageLookup mLineage;

  /**
   * @param lineage array of length equal to the number of samples, giving the the
   * sample id of the predecessor (parent or originating sample) for each sample
   */
  public LineageDenovoChecker(LineageLookup lineage) {
    mLineage = lineage;
  }

  @Override
  public boolean isDenovo(Variant variant, int sample) {
    final VariantSample derived = variant.getSample(sample);
    final int originalId = mLineage.getOriginal(sample);
    final VariantSample original = variant.getSample(originalId);
    final List<String> childAlleles = new ArrayList<>(Arrays.asList(StringUtils.split(derived.getName(), VariantUtils.COLON)));
    for (final String s : StringUtils.split(original.getName(), VariantUtils.COLON)) {
      if (!childAlleles.remove(s)) {
        return true;
      }
    }
    return false;
  }
}
