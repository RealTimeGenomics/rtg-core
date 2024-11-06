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

package com.rtg.segregation;

/**
 * A pattern and whether it is a start of a new phasing region.
 */
final class PatternHolder {

  /**
   * Create a PatternHolder. The input pattern strings indicate which children share which parental haplotype.
   * E.g. <code>001100</code> indicates that children 0,1,4,5 share one haplotype, and 2,3 share the other.
   *
   * @param father a pattern string indicating which children share which paternal haplotype.
   * @param mother a pattern string indicating which children share which maternal haplotype.
   * @param newOrXo a string indicating whether the pattern represents new phase set "N", otherwise a crossover is assumed.
   * @return the PatternHolder
   */
  static PatternHolder fromPatternStrings(String father, String mother, String newOrXo) {
    assert father.length() == mother.length();
    return new PatternHolder(new PatternArray(father, mother), "N".equals(newOrXo));
  }

  private final PatternArray mChP;

  private final boolean mNew;

  PatternHolder(final PatternArray chP, final boolean isnew) {
    mChP = chP;
    mNew = isnew;
  }

  PatternArray pattern() {
    return mChP;
  }

  boolean isNew() {
    return mNew;
  }

  boolean compatible(PatternHolder other) {
    return mChP.compatible(other.mChP);
  }

  String[] incompatibleChildren(PatternHolder other) {
    assert mChP.length() == other.mChP.length();
    assert !compatible(other);
    int minIncompatible = Integer.MAX_VALUE;
    String[] ret = null;
    for (int i = 0; i < Pattern.NUMBER_FLIPS; ++i) {
      int incompatible = 0;
      final String[] compats = new String[mChP.length()];
      for (int j = 0; j < mChP.length(); ++j) {
        if (j >= other.mChP.length()) {
          compats[j] = "N";
        } else if (Pattern.flipIntersect(mChP.index(j), other.mChP.index(j), i) == null) {
          compats[j] = "I";
          ++incompatible;
        } else {
          compats[j] = "C";
        }
      }
      if (incompatible < minIncompatible) {
        minIncompatible = incompatible;
        ret = compats;
      }
    }
    return ret;
  }

  // Get a group number from 0-3 where 0 => 00, 1 => 01, 2 => 10, 3 => 11 from Father:Mother phasing pairs
  int group(int child) {
    int group = 0;
    if ("1".equals(mChP.index(child).faString())) {
      group += 2;
    }
    if ("1".equals(mChP.index(child).moString())) {
      group += 1;
    }
    return group;
  }
}
