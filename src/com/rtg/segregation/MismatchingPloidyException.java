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

import com.rtg.reference.Ploidy;

/**
 * An exception to throw when the Genotype does not match the expected ploidy.
 */
public class MismatchingPloidyException extends Exception {

  /**
   * Constructor
   * @param gt the genotype string
   * @param expectedPloidy the expected ploidy
   */
  public MismatchingPloidyException(String gt, Ploidy expectedPloidy) {
    super("The genotype \"" + gt + "\" does not match the expected ploidy \"" + expectedPloidy + "\"");
  }

}
