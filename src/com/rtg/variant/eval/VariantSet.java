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

package com.rtg.variant.eval;

import java.util.List;
import java.util.Map;

import com.rtg.util.Pair;
import com.rtg.vcf.header.VcfHeader;

/**
 * Iterate over the sequences returning sets of called and baseline variants
 */
public interface VariantSet {

  /**
   * @return the variants for the next sequence or null if there are no more.
   */
  Pair<String, Map<VariantSetType, List<DetectedVariant>>> nextSet();

  /**
   * @return header for baseline files
   */
  VcfHeader baseLineHeader();

  /**
   * @return header for called files
   */
  VcfHeader calledHeader();

  /**
   * @return the number of baseline variants that were skipped during loading
   */
  int getNumberOfSkippedBaselineVariants();

  /**
   * @return the number of called variants that were skipped during loading
   */
  int getNumberOfSkippedCalledVariants();
}
