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
package com.rtg.util;

/**
 * Immutable pair of non-null objects.
 * Care taken with equals and hashcode.
 * @param <A> type of first element in pair.
 * @param <B> type of second element in pair.
 */
public class Pair<A, B> {

  private final A mA;
  private final B mB;

  /**
   * Utility method to make creating pairs easier by skipping typing out the types
   * @param one the first item
   * @param two the second item
   * @param <T1> type of the first item
   * @param <T2> type of the second item
   * @return a pair containing both items
   */
  public static<T1, T2> Pair<T1, T2> create(T1 one, T2 two) {
    return new Pair<>(one, two);
  }

  /**
   * Represents a pair of any type of object.
   * @param a first item
   * @param b second item
   */
  public Pair(final A a, final B b) {
    if (a == null) {
      throw new NullPointerException();
    }
    if (b == null) {
      throw new NullPointerException();
    }
    mA = a;
    mB = b;
  }

  /**
   * Get a.
   * @return Returns the a.
   */
  public A getA() {
    return mA;
  }

  /**
   * Get b.
   * @return Returns the b.
   */
  public B getB() {
    return mB;
  }

  @Override
  public boolean equals(final Object arg0) {
    if (arg0 == null) {
      return false;
    }
    if (arg0 == this) {
      return true;
    }
    if (!(arg0 instanceof Pair)) {
      return false;
    }
    final Pair<?, ?> that = (Pair<?, ?>) arg0;
    return this.mA.equals(that.mA) && this.mB.equals(that.mB);
  }

  @Override
  public int hashCode() {
    return Utils.pairHash(mA.hashCode(), mB.hashCode());
  }

  @Override
  public String toString() {
    return mA + ":" + mB;
  }

}

