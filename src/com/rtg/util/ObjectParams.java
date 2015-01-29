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
package com.rtg.util;

import java.io.IOException;


/**
 * Provides simple (but not incredibly efficient) way of doing hash and equals
 * by using an array of the objects in the class.
 */
public abstract class ObjectParams implements Params {

  /** This is really final but needs to be set late in the constructors for the sub-classes so can't be declared final. */
  protected Object[] mObjects = new Object[0];

  /**
   * Append objects.
   * Useful in subclasses with their own state.
   * @param objs to be added.
   */
  protected void append(final Object[] objs) {
    mObjects = Utils.append(mObjects, objs);
  }

  protected final Object[] objects() {
    return mObjects;
  }

  @Override
  public int hashCode() {
    return Utils.hash(objects());
  }

  @Override
  public boolean equals(final Object obj) {
    if (this == obj) {
      return true;
    }
    if (obj == null) {
      return false;
    }
    assert obj instanceof ObjectParams : obj.getClass().getName();
    final ObjectParams that = (ObjectParams) obj;
    return Utils.equals(this.objects(), that.objects());
  }

  @Override
  public void close() throws IOException {
    // default - do nothing
  }

}
