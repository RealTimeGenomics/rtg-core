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

package com.rtg.vcf.annotation;


/**
 * Annotation for the inbreeding coefficient for a set of samples.
 * @see <a href="http://en.wikipedia.org/wiki/Hardy-Weinberg_principle#Inbreeding_coefficient">Inbreeding coefficient</a>
 */
public class InbreedingCoefficientAnnotation extends AbstractInbreedingCoefficientAnnotation {

  /**
   * Constructor
   */
  public InbreedingCoefficientAnnotation() {
    super("IC", "Inbreeding Coefficient");
  }

  @Override
  protected Double getValue(int total, int hetCount, double expectedHetProbability) {
    if (expectedHetProbability == 0) {
      return null;
    }
    return 1 - hetCount / (total * expectedHetProbability);
  }

}
