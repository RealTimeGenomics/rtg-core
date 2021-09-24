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

package com.rtg.variant.bayes.multisample.population;

import java.util.List;

import com.rtg.variant.bayes.multisample.HaploidDiploidHypotheses;

/**
 */
public interface SiteSpecificPriors {
  /**
   * return new hypotheses if it exists
   * @param templateName template name
   * @param zeroPos position on reference (0 based)
   *
   * @return if population hypotheses exists, a list containing haploid and diploid hypotheses else null
   */
  HaploidDiploidHypotheses<?> getSnpHypotheses(final String templateName, int zeroPos);

  /**
   * Return a list of allele counts for variants in the specified range.
   * @param templateName reference sequence name
   * @param start start of range
   * @param end end of range
   * @return list of <code>AlleleCounts</code>
   */
  List<AlleleCounts> getCounts(final String templateName, int start, int end);

}
