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

import java.util.HashMap;
import java.util.Map;

/**
 * Helps with implementing PseudoEnums.
 * Given the list of values this class will deal with the values and <code>valueOf</code> implementation.
 * @param <T> Class of enum
 */
public class EnumHelper<T extends PseudoEnum> {

  private final Class<T> mClass;

  private final T[] mValues;

  private final String[] mNames;

  private final Map<String, T> mValueOf = new HashMap<>();

  /**
   * Constructs the helper
   * @param c the class of enum this is to help
   * @param values all the enum values in ordinal order.
   */
  public EnumHelper(final Class<T> c, final T[] values) {
    mValues = values.clone();
    mNames = new String[mValues.length];
    int i = 0;
    for (final T t : mValues) {
      mValueOf.put(t.name(), t);
      mNames[i++] = t.name();
    }
    mClass = c;
  }

  /**
   * see {@link java.lang.Enum#valueOf(Class, String)}
   * @param str name of enum
   * @return the enum value
   */
  public T valueOf(final String str) {
    final T ret = mValueOf.get(str);
    if (ret == null) {
      throw new IllegalArgumentException(str + " is not a valid enum value of type: " + mClass.toString());
    }
    return ret;
  }

  /**
   * Like an enum's values() method.
   * @return list of enum values
   */
  public T[] values() {
    return mValues.clone();
  }

  /**
   * Return an array containing the names of all enum members.
   * @return list of enum names
   */
  public String[] names() {
    return mNames.clone();
  }
}
