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

/**
 * Holds a column of data
 */
public abstract class Column implements Cloneable {

  private String mName;

  Column(String name) {
    mName = name;
  }

  public void setName(String name) {
    mName = name;
  }

  /**
   * @return the name of the column
   */
  public String getName() {
    return mName;
  }

  /**
   * @return the number of items
   */
  public abstract int size();

  abstract void add(String value);

  /**
   * Get a value in String representation
   * @param i the item
   * @return the String representation of the item
   */
  abstract String toString(int i);

  abstract Column filter(RegionPredicate p);
}
