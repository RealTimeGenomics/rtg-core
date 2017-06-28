/*
 * Copyright (c) 2017. Real Time Genomics Limited.
 *
 * Use of this source code is bound by the Real Time Genomics Limited Software Licence Agreement
 * for Academic Non-commercial Research Purposes only.
 *
 * If you did not receive a license accompanying this file, a copy must first be obtained by email
 * from support@realtimegenomics.com.  On downloading, using and/or continuing to use this source
 * code you accept the terms of that license agreement and any amendments to those terms that may
 * be made from time to time by Real Time Genomics Limited.
 */

package com.rtg.variant.bayes.complex;

import java.util.List;

import com.rtg.variant.Variant;

/**
 * Declaration of a method for the decomposition of calls into components.
 */
public interface Decomposer {

  /**
   * Given a variant, decompose the given call into smaller constituent parts.
   * An implementation may split the original into an arbitrary number of parts
   * or not split the variant at all.  In addition loci where the original is
   * reference may be trimmed.
   * @param original original call
   * @return decomposed version of the call
   */
  List<Variant> decompose(Variant original);
}
