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
package com.rtg.variant.bayes.multisample.lineage;

import com.rtg.util.Utils;

/**
 * Name of a discrete random variable.
 */
public class Variable implements Comparable<Variable> {

  private final int mSize;
  private final String mName;

  /**
   * Construct a new discrete random variable with a name and size.
   * @param name name of variable
   * @param size number of hypotheses for this variable
   */
  Variable(final String name, final int size) {
    mSize = size;
    mName = name;
  }

  /**
   * Number of entries in the factor.
   * @return size of factor
   */
  int size() {
    return mSize;
  }

  @Override
  public boolean equals(Object that) {
    if (!(that instanceof Variable)) {
      return false;
    }
    final Variable s = (Variable) that;
    return mName.equals(s.mName) && mSize == s.mSize;
  }

  @Override
  public int hashCode() {
    return Utils.pairHash(mName.hashCode(), mSize);
  }

  @Override
  public String toString() {
    return mName;
  }

  @Override
  public int compareTo(Variable o) {
    return mName.compareTo(o.mName);
  }
}
