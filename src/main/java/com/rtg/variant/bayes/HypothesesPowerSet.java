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

package com.rtg.variant.bayes;

import com.rtg.reference.Ploidy;
import com.rtg.variant.bayes.snp.HypothesesNone;
import com.rtg.variant.util.VariantUtils;
import com.rtg.variant.util.arithmetic.PossibilityArithmetic;

/**
 * Provides information about the space of haploid or diploid hypotheses.
 * An hypothesis (whether haploid or diploid) is encoded in a single
 * integer. An hypothesis does not have any associated priors (for that,
 * see <code>HypothesesPrior</code>)
 *
 * @param <D> description type
 */
public class HypothesesPowerSet<D extends Description> extends Hypotheses<D> {

  // Represents a set of hypotheses comprising the power set of the description with
  // the empty set omitted.  Thus, there are 2^{|D|}-1 hypotheses in general and
  // these are denoted internally by the numbers 1, 2, ..., 2^{|D|}-1.  If h is
  // hypothesis, then the bits of h indicate which alleles are present, e.g.
  // h=2^{|D|}-1 has every allele of the description, h=1, h=2, h=4, h=8, etc. has
  // only one allele.  However, for consistency with our other Hypotheses classes
  // externally these hypotheses are numbered 0, ..., 2^{|D|}-2; that means quite
  // a few +/- 1 in this class, but less fiddling outside this class.

  private final int mRef;

  /**
   * Hypotheses which are a power set over the descriptions.
   * @param description of the haploid hypotheses.
   * @param arithmetic to be used for calculations.
   * @param ref the reference allele (in the description)
   */
  public HypothesesPowerSet(final D description, final PossibilityArithmetic arithmetic, final int ref) {
    super(description, arithmetic, new CodePowerSet(description.size()), true);
    mRef = ref;
  }

  // Reference in the hypotheses
  @Override
  public int reference() {
    return mRef == HypothesesNone.NO_HYPOTHESIS ? HypothesesNone.NO_HYPOTHESIS : (1 << mRef) - 1;
  }

  @Override
  public String name(final int hyp) {
    final int code = hyp + 1;
    final StringBuilder sb = new StringBuilder();
    for (int k = 0, v = 1; k < mDescription.size(); ++k, v <<= 1) {
      if ((code & v) != 0) {
        if (sb.length() != 0) {
          sb.append(VariantUtils.COLON);
        }
        sb.append(mDescription.name(k));
      }
    }
    return sb.toString();
  }

  @Override
  public int size() {
    return (1 << mDescription.size()) - 1;
  }

  @Override
  public Ploidy ploidy() {
    return Ploidy.POLYPLOID;
  }

  @Override
  public int maxNameLength() {
    // Longest hypothesis name contain all the alleles
    int len = 0;
    for (int desc = 0; desc < mDescription.size(); ++desc) {
      len += mDescription.name(desc).length();
    }
    return len + mDescription.size() - 1; // add (size - 1) colons
  }
}
