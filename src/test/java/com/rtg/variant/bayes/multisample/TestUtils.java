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

package com.rtg.variant.bayes.multisample;

import com.rtg.variant.Variant;
import com.rtg.variant.VariantLocus;
import com.rtg.variant.VariantSample;

/**
 * Utility methods used in this packages tests.
 */
public final class TestUtils {
  private TestUtils() { }

  /**
   * create a variant give the position
   * @param position variant position
   * @return new variant
   */
  public static Variant createVariant(int position) {
    return createVariant(position, false);
  }

  static Variant createVariant(int position, boolean complex) {
    return createVariant(position, position + 1, complex, 0);
  }

  static Variant createVariant(int position, int end, boolean complex, int indelLength) {
    return createVariant(position, end, complex, indelLength, true);
  }

  static Variant createVariant(int position, int end, boolean complex, int indelLength, boolean interesting) {
    final VariantLocus locus = new VariantLocus("foo", position, end);
    final Variant v = new Variant(locus, new VariantSample[1]);
    if (interesting) {
      v.setInteresting();
    }
    if (position == end || complex) {
      v.setIndel(indelLength);
    }
    return v;
  }

}
