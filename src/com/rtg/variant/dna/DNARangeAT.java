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

package com.rtg.variant.dna;


/**
 */
public final class DNARangeAT {

  private DNARangeAT() { }

  /** Singleton for checking the range. */
  public static final GeneralDNARange DNA = new GeneralDNARange("ACGT", false);

  /** A nucleotide. */
  public static final int A = 0;

  /** C nucleotide. */
  public static final int C = 1;

  /** G nucleotide. */
  public static final int G = 2;

  /** T nucleotide. */
  public static final int T = 3;

}
