/*
 * Copyright (c) $today.theyear. Real Time Genomics Limited.
 *
 * Use of this source code is bound by the Real Time Genomics Limited Software Licence Agreement
 * for Academic Non-commercial Research Purposes only.
 *
 * If you did not receive a license accompanying this file, a copy must first be obtained by email
 * from support@realtimegenomics.com.  On downloading, using and/or continuing to use this source
 * code you accept the terms of that license agreement and any amendments to those terms that may
 * be made from time to time by Real Time Genomics Limited.
 */
package com.rtg.relation;

import java.util.Arrays;

/**
 * Creates a lookup from a sample id to the sample from which it is derived.
 * Somewhat analogous to the ChildFamilyLookup.
 */
public class LineageLookup {

  private final int[] mLineage;

  /**
   * @param lineage the lineage
   */
  public LineageLookup(int... lineage) {
    mLineage = Arrays.copyOf(lineage, lineage.length);
  }

  /**
   * Get the sample id of the sample from which the given sample is derived or -1 if no such orginator exists.
   * @param sampleId id of sample
   * @return original sample id or -1
   */
  public int getOriginal(int sampleId) {
    return mLineage[sampleId];
  }
}
