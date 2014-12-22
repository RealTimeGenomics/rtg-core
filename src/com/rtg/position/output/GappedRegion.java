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
 *         Date: 10/11/11
 *         Time: 3:22 PM
 */
public class GappedRegion extends AbstractGappedRegion<GappedRegion> {
  /**
   *  Construct a concrete implementation of Abstract Gapped region avoiding bizarre type declaration recursion
   */
  GappedRegion(final int id, final BuildParams params, final GapScorer distribution) {
    super(id, params, distribution);
  }
}
