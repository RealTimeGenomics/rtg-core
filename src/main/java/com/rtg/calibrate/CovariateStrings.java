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
  @Override
  public int valueOf(String value) {
    final Integer index = mIndexMap.get(value);
    return index != null ? index : -1;
  }
}
