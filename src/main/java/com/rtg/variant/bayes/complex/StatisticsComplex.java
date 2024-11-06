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
