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
 * A holder for float values - intended to be subclassed as necessary
 * to carry different output formats for rendering etc.
 *
 */
public class FloatValue implements FormattedValue, Comparable<FormattedValue> {

  protected final float mValue;

  FloatValue(final float v) {
    mValue = v;
  }

  @Override
  public String toString() {
    final StringBuilder sb = new StringBuilder();
    toString(sb);
    return sb.toString();
  }


  @Override
  public void toString(final StringBuilder sb) {
    sb.append(mValue);
  }


  @Override
  public int maxLength() {
    return Integer.MAX_VALUE;
  }


  @Override
  public int compareTo(final FormattedValue obj) {
    if (obj == this) {
      return 0;
    }
    if (obj instanceof NullValue) {
      return +1;
    }
    final FloatValue fv = (FloatValue) obj;
    final float v = fv.mValue;
    if (mValue < v) {
      return -1;
    } else if (mValue > v) {
      return +1;
    } else {
      return 0;
    }
  }

  @Override
  public boolean equals(final Object obj) {
    return obj instanceof FloatValue && mValue == ((FloatValue) obj).mValue;
  }

  @Override
  public int hashCode() {
    return Float.floatToRawIntBits(mValue);
  }
}

