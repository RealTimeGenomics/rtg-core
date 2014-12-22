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
package com.rtg.simulation.snpsim;

/**
 * Mutations.
 */
public enum GenDiffMode {
  /** homozygous */
  BOTH_SAME,
  /** heterozygous and both genomes are different */
  DIFFERENT,
  /** heterozygous and only the first genome is different */
  FIRST_ONLY,
  /** heterozygous and only the twin genome is different */
  TWIN_ONLY;

  static String[] diffModeArray(String[] trio, GenDiffMode diffMode) {
    switch (diffMode) {
      case BOTH_SAME:
        return new String[] {trio[0], trio[1]};
      case DIFFERENT:
        return new String[] {trio[0], trio[1], trio[2]};
      case TWIN_ONLY:
        return new String[] {trio[0], trio[0], trio[0].equals(trio[2]) ? trio[1] : trio[2]};
      case FIRST_ONLY:
        return new String[] {trio[0], trio[0].equals(trio[1]) ? trio[2] : trio[1], trio[0]};
      default:
        throw new IllegalStateException("Unpossible");
    }
  }
}
