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
public abstract class Hypotheses<D extends Description> {

  /** Index used to denote invalid hypothesis. */
  public static final int NO_HYPOTHESIS = -1;

  protected final D mDescription;
  private final PossibilityArithmetic mArithmetic;
  private final Code mCode;
  private final boolean mHaploid;

  /**
   * Needed for the cancer case as it does something tricky with {@code code}.
   * @param description of the haploid hypotheses.
   * @param arithmetic to be used for calculations.
   * @param code underlying integer coding.
   * @param haploid true iff the hypotheses are to be haploid.
   */
  protected Hypotheses(final D description, final PossibilityArithmetic arithmetic, final Code code, final boolean haploid) {
    if (description == null || arithmetic == null || code == null) {
      throw new NullPointerException();
    }
    mDescription = description;
    mArithmetic = arithmetic;
    mCode = code;
    mHaploid = haploid;
  }

  /**
   * Common constructor where a Code is automatically selected.
   * @param description of the haploid hypotheses.
   * @param arithmetic to be used for calculations.
   * @param haploid true iff the hypotheses are to be haploid.
   */
  public Hypotheses(final D description, final PossibilityArithmetic arithmetic, final boolean haploid) {
    this(description, arithmetic, haploid ? new CodeHaploid(description.size()) : new CodeDiploid(description.size()), haploid);
  }

  /**
   * @return the description of the haploid hypotheses.
   */
  public D description() {
    return mDescription;
  }

  /**
   * @return the arithmetic used for computing, including priors.
   */
  public PossibilityArithmetic arithmetic() {
    return mArithmetic;
  }

  /**
   * Gets the identity hypothesis.
   * @return the index for the hypothesis that is the same as the reference (if the reference does not correspond to a hypothesis, this should be <code>NOT_A_HYPOTHESIS</code>).
   */
  public abstract int reference();

  /**
   * Gets the name corresponding to one of the hypotheses
   * @param hyp to be displayed.
   * @return human readable string that uniquely describes the hypothesis (usually its nucleotide sequence).
   */
  public String name(int hyp) {
    if (haploid()) {
      return mDescription.name(hyp);
    }
    final int a = code().a(hyp);
    final int b = code().b(hyp);
    //assert b < a;
    return mDescription.name(b) + VariantUtils.COLON + mDescription.name(a);
  }

  /**
   * @return integer coder for hypotheses.
   */
  public Code code() {
    return mCode;
  }

  /**
   * @return the number of valid hypotheses.
   */
  public int size() {
    return code().size();
  }

  /**
   * Check if the index is in valid range.
   * @param index to be checked.
   * @return true iff index is in range.
   */
  public boolean valid(int index) {
    return 0 <= index && index < size();
  }

  /**
   * Determines if the hypotheses are haploid
   * @return true iff the set of hypotheses are haploid.
   */
  public boolean haploid() {
    return mHaploid;
  }

  /**
   * @return the ploidy of the hypotheses.
   */
  public Ploidy ploidy() {
    if (haploid()) {
      return Ploidy.HAPLOID;
    } else {
      return Ploidy.DIPLOID;
    }
  }

  /**
   * Determines if a hypothesis is homozygous
   * @param i hypothesis to be checked.
   * @return true iff i is a homozygous hypothesis.
   */
  public boolean homozygous(final int i) {
    return code().homozygous(i);
  }

  /**
   * @return the maximum length of names.
   */
  public int nameLength() {
    final int haploidLen = mDescription.maxLength();
    final int hypPad;
    if (haploid()) {
      hypPad = haploidLen;
    } else {
      hypPad = 2 * haploidLen + 1;
    }
    return hypPad;
  }


}
