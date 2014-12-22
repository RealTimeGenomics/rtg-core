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

import com.rtg.util.intervals.SequenceNameLocus;
import com.rtg.util.Utils;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;

/**
 * Reference to a variant that selects one side of heterozygous cases and
 * also records a position in the variant chosen.
 */
public class OrientedVariant extends IntegralAbstract implements Comparable<OrientedVariant>, Variant {

  private final Variant mVariant;

  private final boolean mIsAlleleA;

  private double mWeight;

  /**
   * @param variant the variant
   * @param isAlleleA are we taking the A allele
   */
  public OrientedVariant(Variant variant, boolean isAlleleA) {
    super();
    mVariant = variant;
    mIsAlleleA = isAlleleA;
  }

  /**
   * @return the variant.
   */
  Variant variant() {
    return mVariant;
  }


  @Override
  public int compareTo(OrientedVariant that) {
    final int varPos = this.variant().getStart() - that.variant().getStart();
    if (varPos != 0) {
      return varPos;
    }

    if (this.mIsAlleleA == that.mIsAlleleA) {
      return 0;
    }
    return mIsAlleleA ? +1 : -1;
  }

  @Override
  public int hashCode() {
    return Utils.hash(new Object[] {
        mIsAlleleA
        , mVariant
    });
  }

  @Override
  public boolean equals(Object obj) {
    if (obj == null) {
      return false;
    }
    if (obj.getClass() != this.getClass()) {
      return false;
    }
    return this.compareTo((OrientedVariant) obj) == 0;
  }

  @Override
  public void toString(StringBuilder sb) {
    sb.append(mVariant.toString());
    sb.append(" ");
    sb.append(mIsAlleleA ? "+" : "-");
  }

  @Override
  public boolean integrity() {
    Exam.assertNotNull(mVariant);
    Exam.assertNotNull(mVariant.nt(mIsAlleleA));
    return true;
  }

  /**
   * @return true if this is oriented on the A allele
   */
  public boolean isAlleleA() {
    return mIsAlleleA;
  }

  @Override
  public int getStart() {
    return mVariant.getStart();
  }

  @Override
  public int getEnd() {
    return mVariant.getEnd();
  }

  @Override
  public boolean overlaps(SequenceNameLocus other) {
    return mVariant.overlaps(other);
  }

  @Override
  public boolean contains(String sequence, int pos) {
    return mVariant.contains(sequence, pos);
  }

  @Override
  public int getLength() {
    return mVariant.getLength();
  }

  @Override
  public byte[] nt(boolean alleleA) {
    return mVariant.nt(alleleA);
  }

  @Override
  public byte[] ntAlleleA() {
    return mVariant.ntAlleleA();
  }

  @Override
  public byte[] ntAlleleB() {
    return mVariant.ntAlleleB();
  }

  /**
   * @param weight  set the weight
   */
  public void setWeight(double weight) {
    mWeight = weight;
  }

  /**
   * @return calculated weight
   */
  public double getWeight() {
    return mWeight;
  }

  @Override
  public boolean isPhased() {
    return mVariant.isPhased();
  }

  @Override
  public String getSequenceName() {
    return mVariant.getSequenceName();
  }
}
