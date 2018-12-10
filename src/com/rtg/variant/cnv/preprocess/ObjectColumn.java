/*
 * Copyright (c) 2016. Real Time Genomics Limited.
 *
 * Use of this source code is bound by the Real Time Genomics Limited Software Licence Agreement
 * for Academic Non-commercial Research Purposes only.
 *
 * If you did not receive a license accompanying this file, a copy must first be obtained by email
 * from support@realtimegenomics.com.  On downloading, using and/or continuing to use this source
 * code you accept the terms of that license agreement and any amendments to those terms that may
 * be made from time to time by Real Time Genomics Limited.
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
