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
package com.rtg.simulation;

/**
 * Class stores the results for mappings
 */
class ReadSimEvalStatistics {

  private static final byte FOUND = 0x01;
  private static final byte MATED = 0x02;
  private static final byte UNMATED = 0x04;
  private static final byte UNMAPPED = 0x08;
  private static final byte MULTIPLE = 0x10;
  private static final byte MAPPED = 0x20;

  private final int mLen;
  private final byte[] mData;

  ReadSimEvalStatistics(int i) {
    mLen = i;
    mData = new byte[mLen];
  }

  int length() {
    return mLen;
  }

  void found(int i) {
    set(i, FOUND);
  }

  boolean isFound(int i) {
    return get(i, FOUND);
  }

  void mated(int i) {
    set(i, MATED);
  }

  boolean isMated(int i) {
    return get(i, MATED);
  }

  void unmated(int i) {
    set(i, UNMATED);
  }

  boolean isUnmated(int i) {
    return get(i, UNMATED);
  }

  void unmapped(int i) {
    set(i, UNMAPPED);
  }

  boolean isUnmapped(int i) {
    return get(i, UNMAPPED);
  }

  void mapped(int i) {
    set(i, MAPPED);
  }

  boolean isMapped(int i) {
    return get(i, MAPPED);
  }

  void multiple(int i) {
    set(i, MULTIPLE);
  }

  boolean isMultiple(int i) {
    return get(i, MULTIPLE);
  }

  private boolean get(int readId, byte flag) {
    return (mData[readId] & flag) == flag;
  }

  private void set(int readId, byte flag) {
    mData[readId] |= flag;
  }

}
