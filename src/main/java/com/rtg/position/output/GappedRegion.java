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
package com.rtg.position.output;

import com.rtg.launcher.BuildParams;

/**
 */
public class GappedRegion extends AbstractGappedRegion<GappedRegion> {
  /**
   * Construct a concrete implementation of Abstract Gapped region avoiding bizarre type declaration recursion
   * @param id an identifier unique within a particular <code>GappedOutput</code> object.
   * @param params parameters for the query.
   * @param distribution probability distribution for gap probabilities.
   */
  GappedRegion(final int id, final BuildParams params, final GapScorer distribution) {
    super(id, params, distribution);
  }
}
