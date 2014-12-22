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
package com.rtg.util.format;

/**
 * A holder for null or uncomputable values.
 *
 */
public final class NullValue implements FormattedValue, Comparable<FormattedValue> {

  /** Singleton null value. */
  public static final NullValue NULL = new NullValue();

  //private so only instace is the static above
  private NullValue() { }

  @Override
  public String toString() {
    return "--";
  }


  @Override
  public void toString(final StringBuilder sb) {
    sb.append("--");
  }


  @Override
  public int maxLength() {
    return 2;
  }


  @Override
  public int compareTo(final FormattedValue obj) {
    if (obj == this) {
      return 0;
    }
    return -1; //it is less than everything else;
  }

  @Override
  public boolean equals(final Object obj) {
    return obj == this;
  }

  @Override
  public int hashCode() {
    return -1;
  }
}

