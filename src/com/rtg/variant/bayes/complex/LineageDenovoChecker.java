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

import com.rtg.util.StringUtils;
import com.rtg.variant.Variant;
import com.rtg.variant.VariantSample;
import com.rtg.variant.bayes.multisample.lineage.Lineage;

/**
 * De novo checker that supports cell lineage de novo checking
 */
public class LineageDenovoChecker implements DenovoChecker {

  final Lineage mLineage;

  /**
   * @param lineage the cell lineage
   */
  public LineageDenovoChecker(Lineage lineage) {
    mLineage = lineage;
  }

  @Override
  public boolean isDenovo(Variant variant, int sample) {
    final VariantSample child = variant.getSample(sample);
    final int parentId = mLineage.parent(sample);
    final VariantSample parent = variant.getSample(parentId);
    final List<String> childAlleles = new ArrayList<>(Arrays.asList(StringUtils.split(child.getName(), ':')));
    final List<String> parentAlleles = new ArrayList<>(Arrays.asList(StringUtils.split(parent.getName(), ':')));
    for (String s : parentAlleles) {
      if (!childAlleles.remove(s)) {
        return true;
      }
    }
    return false;
  }
}
