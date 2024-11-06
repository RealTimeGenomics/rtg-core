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
package com.rtg.variant.cnv.preprocess;

import java.util.ArrayList;
import java.util.Collection;

import com.reeltwo.jumble.annotations.TestClass;

/**
 * Holds a column of data.
 */
@TestClass("com.rtg.variant.cnv.preprocess.StringColumnTest")
public abstract class ObjectColumn<T> extends Column {

  private ArrayList<T> mData = new ArrayList<>();

  ObjectColumn(String name) {
    super(name);
  }

  @Override
  public int size() {
    return mData.size();
  }

  @Override
  String toString(int i) {
    return String.valueOf(get(i));
  }

  /**
   * @return the values
   */
  public Iterable<T> values() {
    return mData;
  }

  /**
   * Add a new entry
   * @param value the value to add
   */
  void add(T value) {
    mData.add(value);
  }

  /**
   * Add a value at a specific index
   * @param i the index of the entry to insert at
   * @param value the value to add
   */
  void add(int i, T value) {
    mData.add(i, value);
  }

  void set(Collection<T> values) {
    mData = new ArrayList<>();
    mData.addAll(values);
  }

  /**
   * Set the value of an entry
   * @param i the index of the entry to set
   * @param value the value to add
   */
  void set(int i, T value) {
    mData.set(i, value);
  }

  /**
   * @param i the index of the value to get
   * @return the value of an entry
   */
  public T get(int i) {
    return mData.get(i);
  }

  @Override
  public void remove(int i) {
    mData.remove(i);
  }

  @Override
  @SuppressWarnings("unchecked")
  ObjectColumn<T> filter(RegionPredicate p) {
    final ArrayList<T> newData = new ArrayList<>();
    for (int i = 0; i < size(); ++i) {
      if (p.test(i)) {
        newData.add(mData.get(i));
      }
    }
    try {
      final ObjectColumn<T> result = (ObjectColumn<T>) clone();
      result.set(newData);
      return result;
    } catch (CloneNotSupportedException e) {
      throw new RuntimeException(e);
    }
  }

}
