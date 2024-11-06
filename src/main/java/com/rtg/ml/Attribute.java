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
package com.rtg.ml;

import java.util.ArrayList;
import java.util.HashMap;

/**
 * Encapsulate an attribute.
 *
 */
public final class Attribute {

  private final String mName;
  private final MlDataType mDataType;
  private final ValueEncoder mEncoder;

  /**
   * Construct a container for a dataset with specified attributes.
   * @param name name of the attribute
   * @param mlDataType data type for attribute
   */
  public Attribute(final String name, MlDataType mlDataType) {
    if (name == null) {
      throw new NullPointerException("Null name");
    }
    mName = name;
    mDataType = mlDataType;
    switch (mDataType) {
      case STRING:
        mEncoder = new StringEncoder();
        break;
      case BOOLEAN:
        mEncoder = new BooleanEncoder();
        break;
      case INTEGER:
        mEncoder = new IntEncoder();
        break;
      case DOUBLE:
        mEncoder = new DoubleEncoder();
        break;
      default:
        throw new UnsupportedOperationException();
    }
  }

  public String getName() {
    return mName;
  }

  public MlDataType getDataType() {
    return mDataType;
  }

  /**
   * @param encodedValue value encoded using {@code encodeValue} method
   * @return true if encoded value is the missing value (null)
   */
  public static boolean isMissingValue(double encodedValue) {
    return Double.isNaN(encodedValue);
  }

  /**
   * Encode the object value into a consistent data type for use within classifiers
   * @param val the original value
   * @return the encoded value
   */
  public double encodeValue(Object val) {
    if (val == null) {
      return Double.NaN;
    }
    return mEncoder.encode(val);
  }

  /**
   * Decode the value back into its original form
   * @param val the encoded value
   * @return the original value
   */
  public Object decodeValue(double val) {
    if (isMissingValue(val)) {
      return null;
    }
    return mEncoder.decode(val);
  }

  /**
   * @return number of distinct values possible/seen in a nominal value (string or boolean type)
   */
  public int nominalSize() {
    if (mDataType.isNumeric()) {
      throw new UnsupportedOperationException("Enumeration only possible on nominal data types");
    }
    return mEncoder.nominalSize();
  }

  /**
   * Return an ordered comparison over the underlying values. This may be a relatively expensive
   * operation if internal decoding is required.
   * @param a1 first value to compare
   * @param a2 second value to compare
   * @return integer indicating ordering
   */
  public int compare(double a1, double a2) {
    if (a1 == a2) {
      return 0;
    }
    if (Attribute.isMissingValue(a1)) {
      return -1;
    } else if (Attribute.isMissingValue(a2)) {
      return 1;
    }
    return mEncoder.compare(a1, a2);
  }

  private interface ValueEncoder {
    double encode(Object value);
    Object decode(double value);
    int nominalSize();
    int compare(double value1, double value2);
  }
  private static final class IntEncoder implements ValueEncoder {
    @Override
    public double encode(Object value) {

      return (Integer) value;
    }

    @Override
    public Object decode(double value) {
      return (int) value;
    }

    @Override
    public int compare(double a1, double a2) {
      return Double.compare(a1, a2);
    }

    @Override
    public int nominalSize() {
      throw new UnsupportedOperationException("Not implemented yet");
    }
  }
  private static final class DoubleEncoder implements ValueEncoder {
    @Override
    public double encode(Object value) {
      return (Double) value;
    }

    @Override
    public Object decode(double value) {
      return value;
    }

    @Override
    public int compare(double a1, double a2) {
      return Double.compare(a1, a2);
    }

    @Override
    public int nominalSize() {
      throw new UnsupportedOperationException("Not implemented yet");
    }
  }
  private static final class StringEncoder implements ValueEncoder {
    // Note that a particular mapping of String value -> encoded value is transient and depends on the order in which
    // values are encoded.
    // In serialized models we store unencoded values and the mapping is repopulated during deserialization.
    private final HashMap<String, Integer> mValues = new HashMap<>();
    private final ArrayList<String> mDecoder = new ArrayList<>();

    @Override
    public synchronized double encode(Object value) {
      final String strVal = (String) value;
      Integer encoded = mValues.get(strVal);
      if (encoded == null) {
        encoded = mDecoder.size();
        mValues.put(strVal, encoded);
        mDecoder.add(strVal);
      }
      return (double) encoded;
    }

    @Override
    public synchronized Object decode(double value) {
      final int index = (int) value;
      return mDecoder.get(index);
    }

    @Override
    public synchronized int nominalSize() {
      return mDecoder.size();
    }

    @Override
    @SuppressWarnings(value = {"unchecked", "rawtypes"})
    public int compare(double a1, double a2) {
      return ((Comparable) decode(a1)).compareTo(decode(a2));
    }
  }
  private static final class BooleanEncoder implements ValueEncoder {
    @Override
    public double encode(Object value) {
      return ((Boolean) value) ? 1.0 : 0.0;
    }

    @Override
    public Object decode(double value) {
      return value > 0.0;
    }

    @Override
    public int compare(double a1, double a2) {
      return Double.compare(a1, a2);
    }

    @Override
    public int nominalSize() {
      return 2;
    }
  }
}
