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

package com.rtg.variant;

import com.rtg.sam.BadSuperCigarException;

/**
 */
public interface ReadParserInterface {

  /**
   * Expands read input and calls appropriate matcher.
   * @param me machine errors used for estimating error rates in calling.
   * @param var variant alignment record to process
   * @param qdefault if quality is null default quality value to use.
   * @param templateBytes template array with nucleotide codes (coded as 0 = N ... 4 = T).
   * @throws BadSuperCigarException if a recognized problem with a cigar is found.
   */
  void toMatcher(AbstractMachineErrorParams me, VariantAlignmentRecord var, int qdefault, byte[] templateBytes) throws BadSuperCigarException;

}
