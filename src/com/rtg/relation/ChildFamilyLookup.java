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
package com.rtg.relation;

/**
 * Creates a lookup from a sample id to the family in which that sample is a child
 */
public class ChildFamilyLookup {
  private final Family[] mChildToFamily;

  /**
   * @param numSamples total number of samples
   * @param families the set of families, with sample ids set appropriately
   */
  public ChildFamilyLookup(int numSamples, Family... families) {
    final Family[] childToFamily = new Family[numSamples];
    if (families != null) {
      for (Family f : families) {
        final int[] sampleIds = f.getSampleIds();
        for (int i = Family.FIRST_CHILD_INDEX; i < sampleIds.length; i++) {
          childToFamily[sampleIds[i]] = f;
        }
      }
    }
    mChildToFamily = childToFamily;
  }

  /**
   * Get the family (if any) in which given sample is a child.
   * @param sampleId id of sample
   * @return the family in which the sample is a child, or null if given sample is missing one or both parents from the genome relationships input
   */
  public Family getFamily(int sampleId) {
    return mChildToFamily[sampleId];
  }
}
