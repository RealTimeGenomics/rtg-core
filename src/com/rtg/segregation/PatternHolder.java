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
