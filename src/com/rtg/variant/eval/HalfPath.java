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

package com.rtg.variant.eval;

import com.rtg.util.BasicLinkedListNode;
import com.rtg.util.Utils;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;

/**
 * One half of a path that reconciles two sequences of variants.
 */
public class HalfPath extends IntegralAbstract implements Comparable<HalfPath> {

  private final HaplotypePlayback mHaplotypeA;
  private HaplotypePlayback mHaplotypeB;

  private BasicLinkedListNode<OrientedVariant> mIncluded;
  private BasicLinkedListNode<Variant> mExcluded;

  private int mVariantEndPosition;
  Variant mLastVariant = null;

  private boolean mFinishedTypeA;
  private boolean mFinishedTypeB;

  /**
   * Construct an empty <code>HalfPath</code>
   * @param template the template on which the half path resides
   */
  HalfPath(byte[] template) {
    mHaplotypeA = new HaplotypePlayback(template);
    mHaplotypeB = null;
  }

  /**
   * Construct a child <code>HalfPath</code>
   * @param path the parent half path
   */
  HalfPath(HalfPath path) {
    mIncluded = path.mIncluded;
    mExcluded = path.mExcluded;
    mHaplotypeA = path.mHaplotypeA.copy();
    mHaplotypeB = path.mHaplotypeB == null ? null : path.mHaplotypeB.copy();
    mVariantEndPosition = path.mVariantEndPosition;
    mLastVariant = path.mLastVariant;
    mFinishedTypeA = path.mFinishedTypeA;
    mFinishedTypeB = path.mFinishedTypeB;
  }

  /**
   * @param var the variant to check
   * @return true if the variant is later than the last included one.
   */
  public boolean isNew(Variant var) {
    return mLastVariant == null || mLastVariant.getStart() < var.getStart() || mLastVariant.getEnd() < var.getEnd();
  }

  void exclude(Variant var) {
    mExcluded = new BasicLinkedListNode<>(var, mExcluded);
    mVariantEndPosition = var.getEnd();
    mLastVariant = var;
  }

  void include(OrientedVariant var) {
    mIncluded = new BasicLinkedListNode<>(var, mIncluded);
    mVariantEndPosition = var.variant().getEnd();
    mLastVariant = var.variant();

    // Lazy creation of the B haplotype
    if (mHaplotypeB == null && var.variant().ntAlleleB() != null) {
      mHaplotypeB = mHaplotypeA.copy();
    }

    mHaplotypeA.addVariant(var);

    if (mHaplotypeB != null) {
      if (var.variant().ntAlleleB() == null) { // Homozygous variant
        mHaplotypeB.addVariant(var);
      } else { // Add the other side of the heterozygous variant.
        mHaplotypeB.addVariant(new OrientedVariant(var.variant(), !var.isAlleleA()));
      }
    }

  }

  void step() {
    haplotypeAStep();
    haplotypeBStep();
  }

  void haplotypeAStep() {
    if (mHaplotypeB == null) {
      mHaplotypeB = mHaplotypeA.copy();
    }
    if (mHaplotypeA.hasNext()) {
      mHaplotypeA.next();
    } else {
      mFinishedTypeA = true;
    }
  }

  void haplotypeBStep() {
    if (mHaplotypeB == null) {
      mHaplotypeB = mHaplotypeA.copy();
    }
    if (mHaplotypeB.hasNext()) {
      mHaplotypeB.next();
    } else {
      mFinishedTypeB = true;
    }
  }

  /**
   * @return true if the half path has reached the end of the template on both haplotypes
   */
  public boolean finished() {
    return mFinishedTypeB && mFinishedTypeA;
  }

  /**
   * @return true if the half path has reached the end of the template on the A haplotype
   */
  public boolean finishedHaplotypeA() {
    return mFinishedTypeA;
  }

  /**
   * @return true if the half path has reached the end of the template on the B haplotype
   */
  public boolean finishedHaplotypeB() {
    return mFinishedTypeB;
  }

  public BasicLinkedListNode<OrientedVariant> getIncluded() {
    return mIncluded;
  }

  public BasicLinkedListNode<Variant> getExcluded() {
    return mExcluded;
  }

  /**
   * Get variant end position.
   * @return Returns the variant end position.
   */
  public int getVariantEndPosition() {
    return mVariantEndPosition;
  }

  /**
   * Retrieve the leading reference so we can add the next variants
   * @return the leading reference
   */
  public int getPosition() {
    if (mHaplotypeB == null) {
      return mHaplotypeA.templatePosition();
    }
    if (mHaplotypeA.templatePosition() > mHaplotypeB.templatePosition()) {
      return mHaplotypeA.templatePosition();
    } else {
      return mHaplotypeB.templatePosition();
    }
  }

  /**
   * @return the next base on the A haplotype
   */
  public byte nextHaplotypeABase() {
    return nextBase(true);
  }

  /**
   * @return the next base on the B haplotype
   */
  public byte nextHaplotypeBBase() {
    return nextBase(false);
  }

  private byte nextBase(boolean haplotypeA) {
    if (haplotypeA || mHaplotypeB == null) {
      return mHaplotypeA.nt();
    } else {
      return mHaplotypeB.nt();
    }
  }
  void moveForward(int position) {
    if (mHaplotypeB != null) {
      mHaplotypeB.moveForward(position);
    }
    mHaplotypeA.moveForward(position);
  }

  @Override
  public int compareTo(HalfPath that) {
    final int c0 = this.mHaplotypeA.templatePosition() - that.mHaplotypeA.templatePosition();
    if (c0 != 0) {
      return c0;
    }

    if (this.mHaplotypeB == null && that.mHaplotypeB != null) {
      return -1;
    }
    if (this.mHaplotypeB != null && that.mHaplotypeB == null) {
      return +1;
    }
    final int plus = this.mHaplotypeA.compareTo(that.mHaplotypeA);
    if (this.mHaplotypeB == null || plus != 0) {
      return plus;
    }
    return this.mHaplotypeB.compareTo(that.mHaplotypeB);
  }

  /**
   * performs a compare of the template position of both haplotypes
   * @return &lt; 0 if the B haplotype is leading, &gt; 0  if A haplotype is leading, 0 if in the same place
   */
  public int compareHaplotypePositions() {
    if (mHaplotypeB == null) {
      return 0;
    }
    return mHaplotypeA.templatePosition() - mHaplotypeB.templatePosition();
  }


  @Override
  public boolean equals(Object obj) {
    return obj != null && compareTo((HalfPath) obj) == 0;
  }

  @Override
  public void toString(StringBuilder sb) {
    sb.append(" +:");
    sb.append(mHaplotypeA.templatePosition());
    if (mHaplotypeB != null) {
      sb.append(" -:");
      sb.append(mHaplotypeB.templatePosition());
    }
    sb.append(LS);
    sb.append("included:");
    sb.append(listToString(mIncluded));
    sb.append(LS);
    sb.append("excluded:");
    sb.append(listToString(mExcluded));
    sb.append(LS);
  }
  @Override
  public int hashCode() {
    return Utils.hash(new Object[] {
      mHaplotypeB
      , mHaplotypeA
      , mVariantEndPosition
    });
  }

  static <T> String listToString(BasicLinkedListNode<T> list) {
    final StringBuilder sb = new StringBuilder();
    sb.append("[");
    String join = "";
    if (list != null) {
      for (final T v : list) {
        sb.append(join);
        sb.append(v.toString());
        join = ", ";
      }
    }
    sb.append("]");
    return sb.toString();
  }

  @Override
  public boolean integrity() {
    Exam.assertNotNull(mHaplotypeA);

    return true;
  }

}
