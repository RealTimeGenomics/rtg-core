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

import junit.framework.TestCase;

/**
 */
public class MismatchingPloidyExceptionTest extends TestCase {

  public void test() {
    final MismatchingPloidyException ex = new MismatchingPloidyException("0/1", Ploidy.HAPLOID);
    assertEquals("The genotype \"0/1\" does not match the expected ploidy \"HAPLOID\"", ex.getMessage());
  }
}
