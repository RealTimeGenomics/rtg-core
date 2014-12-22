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

package com.rtg.variant;

import com.rtg.relation.Relationship;
import com.rtg.relation.Relationship.RelationshipType;
import com.rtg.util.StringUtils;

/**
 */
public class SomaticStatistics extends VariantStatistics {

  private final String mDerivedName;

  /**
   * Somatic specific variant statistics which only outputs the statistics for the derived sample.
   * @param params the params
   */
  public SomaticStatistics(VariantParams params) {
    super(params.directory());
    final Relationship[] derived = params.genomeRelationships().relationships(RelationshipType.ORIGINAL_DERIVED);
    final int numberOfGenomes = params.genomeRelationships().genomes().length;
    assert derived.length == 1 || numberOfGenomes == 2;
    final String[] genomeNames = {derived[0].first(), derived[0].second()};
    mDerivedName = genomeNames[1] ;
  }

  @Override
  public String getStatistics() {
    final StringBuilder out = new StringBuilder();
    out.append("Passed Filters               : ").append(mTotalVariants).append(StringUtils.LS);
    out.append("Failed Filters               : ").append(mTotalFiltered).append(StringUtils.LS);
    if (mExcessCoverage > 0) {
      out.append("Excessive Coverage           : ").append(mExcessCoverage).append(StringUtils.LS);
    }
    if (mExcessHypotheses > 0) {
      out.append("Excessive Hypotheses         : ").append(mExcessHypotheses).append(StringUtils.LS);
    }
    if (mNoHypotheses > 0) {
      out.append("No Hypotheses                : ").append(mNoHypotheses).append(StringUtils.LS);
    }

    final PerSampleVariantStatistics derivedSampleStats = mPerSampleStats.get(mDerivedName);
    if (derivedSampleStats != null) {
      derivedSampleStats.appendStatistics(out);
    }

    return out.toString();
  }
}
