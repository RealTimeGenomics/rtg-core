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

package com.rtg.calibrate;

import com.reeltwo.jumble.annotations.TestClass;

/**
 */
@TestClass({"com.rtg.calibrate.CovariateReadGroupTest", "com.rtg.calibrate.CovariateBaseQualityTest"})
public abstract class CovariateImpl implements Covariate {

  private final String mName;

  private int mSize;
  protected int mNewSize;

  /**
   * Default constructor for <code>CovariateImpl</code>
   * @param name the name of the variable (lower case)
   * @param size the initial size of the variable
   */
  public CovariateImpl(String name, int size) {
    assert size >= 0;
    mSize = mNewSize = size;
    mName = name;
  }

  @Override
  public String name() {
    return mName;
  }

  @Override
  public final int size() {
    return mSize;
  }

  @Override
  public void resized() {
    mSize = newSize();
  }

  @Override
  public final int newSize() {
    return mNewSize;
  }

  protected final void setNewSize(int size) {
    assert size >= mNewSize;
    mNewSize = size;
  }

  @Override
  public final boolean sizeChanged() {
    assert newSize() >= mSize;
    return newSize() > mSize;
  }

  /**
   * Default implementation returns number as string.
   */
  @Override
  public String valueString(int value) {
    return Integer.toString(value);
  }

  /**
   * Default implementation returns string parsed as number.
   */
  @Override
  public int parse(String value) {
    final int ret = Integer.parseInt(value);
    if (ret >= newSize()) {
      setNewSize(ret + 1);
    }
    return ret;
  }

  /**
   * Default implementation returns string parsed as number.
   */
  @Override
  public int valueOf(String value) {
    return parse(value);
  }

  @Override
  public String toString() {
    return name() + ":" + size() + ":" + newSize();
  }


}
