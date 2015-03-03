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
package com.rtg.util.memory;

import java.util.Stack;

/**
 * An <code>ObjectFilter</code> that rejects instances of particular
 * classes. This implementation only compares exact classes, not
 * interfaces or subclasses.
 */
public class ExcludeClassFilter {

  private IdentitySet mExcludeSet = null;

  /**
   * Creates a new <code>ExcludeClassFilter</code> instance.
   *
   * @param excludeClasses classes to exclude
   */
  public ExcludeClassFilter(final Class<?>[] excludeClasses) {
    setExcludedClasses(excludeClasses);
  }

  /**
   * Describe <code>setExcludedClasses</code> method here.
   *
   * @param excludeClasses classes to exclude
   */
  public void setExcludedClasses(final Class<?>[] excludeClasses) {
    if (excludeClasses == null) {
      mExcludeSet = null;
    } else {
      mExcludeSet = new IdentitySet();
      mExcludeSet.addAll(excludeClasses);
    }
  }

  /**
   * Returns true if the object is acceptable.
   *
   * @param path a <code>Stack</code> containing the path of objects
   * walked to reach this object.
   * @param obj an <code>Object</code> to test
   * @return true if the object is accepted.
   */
  public boolean accept(final Stack<Object> path, final Object obj) {
    return mExcludeSet == null || !mExcludeSet.contains(obj.getClass());
  }
}

