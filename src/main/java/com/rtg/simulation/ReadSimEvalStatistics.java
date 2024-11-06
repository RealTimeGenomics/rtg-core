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
