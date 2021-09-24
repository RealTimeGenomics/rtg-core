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

import com.rtg.variant.bayes.Description;
import com.rtg.variant.bayes.EvidenceInterface;
import com.rtg.variant.bayes.StatisticsDouble;
import com.rtg.variant.match.Match;

/**
 * Statistics specific to the complex caller.
 */
public class StatisticsComplex extends StatisticsDouble {

  private final int mRefLength;

  /**
   * @param description the description
   * @param refLength the reference length (used for coverage calculation for partially overlapping evidence)
   */
  public StatisticsComplex(Description description, final int refLength) {
    super(description);
    mRefLength = refLength;
  }

  @Override
  protected double coverageIncrement(EvidenceInterface distribution) {
    if (distribution instanceof EvidenceComplex) {
      return coverageIncrement(((EvidenceComplex) distribution).match());
    }
    return super.coverageIncrement(distribution);
  }

  @Override
  protected double errorIncrement(EvidenceInterface evidence) {
    if (evidence instanceof EvidenceComplex) {
      return ((EvidenceComplex) evidence).match().correction() * coverageIncrement(((EvidenceComplex) evidence).match());
    }
    return super.errorIncrement(evidence);
  }


  double coverageIncrement(Match m) {
    final double ret;
    if (mRefLength == 0 || m.isFixedLeft() && m.isFixedRight()) {
      ret = 1.0;
    } else {
      final double b = (double) m.length() / mRefLength;
      ret = Math.min(1.0, b);
    }
    return ret;
  }

}
