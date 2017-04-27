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
package com.rtg.ml;

import java.util.ArrayList;
import java.util.HashMap;

/**
 * Encapsulate an attribute.
 *
 */
public class Attribute {

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


  private interface ValueEncoder {
    double encode(Object value);
    Object decode(double value);
    int nominalSize();
  }
  private static class IntEncoder implements ValueEncoder {
    @Override
    public double encode(Object value) {

      return (Integer) value;
    }

    @Override
    public Object decode(double value) {
      return (int) value;
    }

    @Override
    public int nominalSize() {
      throw new UnsupportedOperationException("Not implemented yet");
    }
  }
  private static class DoubleEncoder implements ValueEncoder {
    @Override
    public double encode(Object value) {
      return (Double) value;
    }

    @Override
    public Object decode(double value) {
      return value;
    }

    @Override
    public int nominalSize() {
      throw new UnsupportedOperationException("Not implemented yet");
    }
  }
  private static class StringEncoder implements ValueEncoder {
    // Note that a particular mapping of String value -> encoded value is transient and depends on the order in which
    // values are encoded.
    // In serialized models we store unencoded values and the mapping is repopulated during deserialization.
    private final HashMap<String, Integer> mValues = new HashMap<>();
    private final ArrayList<String> mDecoder = new ArrayList<>();

    @Override
    public double encode(Object value) {
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
    public Object decode(double value) {
      final int index = (int) value;
      return mDecoder.get(index);
    }

    @Override
    public int nominalSize() {
      return mDecoder.size();
    }
  }
  private static class BooleanEncoder implements ValueEncoder {
    @Override
    public double encode(Object value) {
      return ((Boolean) value) ? 1.0 : 0.0;
    }

    @Override
    public Object decode(double value) {
      return value > 0.0;
    }

    @Override
    public int nominalSize() {
      return 2;
    }
  }
}
