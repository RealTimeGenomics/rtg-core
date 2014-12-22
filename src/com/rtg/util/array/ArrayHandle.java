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
package com.rtg.util.array;

/**
 * Holds an <code>ArrayType</code> and a length.
 * Can be used to create an array.
 */
public class ArrayHandle {

  private final ArrayType mType;

  private final long mLength;

  /**
   * Construct a handle.
   * @param type of underlying index.
   * @param length of index.
   */
  public ArrayHandle(final ArrayType type, final long length) {
    if (type == null || length < 0) {
      throw new IllegalArgumentException();
    }
    mType = type;
    mLength = length;
  }

  /**
   * Create  an unsigned index of the specified length and type.
   * @return an unsigned index of the specified length and type.
   */
  public ExtensibleIndex createUnsigned() {
    return mType.createUnsigned(mLength);
  }

  /**
   * Get  the length.
   * @return the length.
   */
  public long length() {
    return mLength;
  }

  /**
   * Get the type.
   * @return the type.
   */
  public ArrayType type() {
    return mType;
  }

  /**
   * Compute the number of bytes consumed by the created arrays.
   * @return the number of bytes consumed by the created arrays.
   */
  public long bytes() {
    return mType.bytes(mLength);
  }

  @Override
  public String toString() {
    return "";
  }

}

