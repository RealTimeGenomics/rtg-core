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

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.reeltwo.jumble.annotations.TestClass;

/**
 * Base implementation of a String based variable
 */
@TestClass({"com.rtg.calibrate.CovariateReadGroupTest", "com.rtg.calibrate.CovariateSequenceTest"})
public abstract class CovariateStrings extends CovariateImpl {

  private final List<String> mStrings = new ArrayList<>();
  private final Map<String, Integer> mIndexMap = new HashMap<>();

  /**
   * Default constructor for <code>CovariateStrings</code>
   * @param name the name of the variable (lower case)
   */
  public CovariateStrings(String name) {
    super(name, 0);
  }

  @Override
  public String valueString(int value) {
    return mStrings.get(value);
  }

  private int add(String value) {
    mStrings.add(value);
    final int size = newSize();
    final int listsize = mStrings.size();
    if (listsize > size) {
      setNewSize(size + selectIncrement());
    }
    final int index = listsize - 1;
    mIndexMap.put(value, index);
    return index;
  }

  /**
   * Unknown strings are added automatically.
   */
  @Override
  public int parse(String value) {
    final Integer pos = mIndexMap.get(value);
    if (pos != null) {
      return pos;
    } else {
      return add(value);
    }
  }

  private int selectIncrement() {
    final int size = newSize();
    if (size >= 1000) {
      return 1000;
    } else if (size >= 100) {
      return 100;
    } else if (size >= 10) {
      return 10;
    }
    return 1;
  }

  /**
   */
  @Override
  public int valueOf(String value) {
    final Integer index = mIndexMap.get(value);
    return index != null ? index : -1;
  }
}
