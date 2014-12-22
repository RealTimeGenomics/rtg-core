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

import com.rtg.util.intervals.SequenceNameLocus;

/**
 * Hold a variant.
 */
public interface Variant extends SequenceNameLocus {

  /**
   * One allele of the variant as determined by allele parameter.
   * @param alleleA if true select the A allele.
   * @return the allele (may be null or zero length)
   */
  byte[] nt(boolean alleleA);

  /**
   * The first allele of the variant. Always non null but may be zero length.
   * @return the A allele
   */
  byte[] ntAlleleA();

  /**
   * Second allele of the variant. Will be null if homozygous and may be zero length.
   * @return the B allele, or null if homozygous
   */
  byte[] ntAlleleB();

  /**
   * @return true if the call has phasing information
   */
  boolean isPhased();
}
