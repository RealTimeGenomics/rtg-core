/*
 * Copyright (c) 2015. Real Time Genomics Limited.
 *
 * Use of this source code is bound by the Real Time Genomics Limited Software Licence Agreement
 * for Academic Non-commercial Research Purposes only.
 *
 * If you did not receive a license accompanying this file, a copy must first be obtained by email
 * from support@realtimegenomics.com.  On downloading, using and/or continuing to use this source
 * code you accept the terms of that license agreement and any amendments to those terms that may
 * be made from time to time by Real Time Genomics Limited.
 */
package com.rtg.variant.eval;

import java.util.List;

import com.reeltwo.jumble.annotations.JumbleIgnore;

/**
 * A simple holder for loaded results plus stats about how many variants were skipped.
 */
@JumbleIgnore
class LoadedVariants {

  final int mSkippedDuringLoading;
  final List<DetectedVariant> mVariants;

  /**
   * Construct the holder
   * @param variants the variants that were loaded
   * @param skipped a count of the number of variants skipped during loading
   */
  public LoadedVariants(List<DetectedVariant> variants, int skipped) {
    mVariants = variants;
    mSkippedDuringLoading = skipped;
  }
}
