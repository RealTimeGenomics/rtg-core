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

package com.rtg.variant.bayes;

import com.rtg.reference.Ploidy;
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
    return (1 << mRef) - 1;
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
