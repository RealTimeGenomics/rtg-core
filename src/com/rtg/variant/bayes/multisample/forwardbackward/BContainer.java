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
package com.rtg.variant.bayes.multisample.forwardbackward;

import com.rtg.variant.bayes.Factor;

/**
 * Contains B values for an individual sample (for definition of B see <code>pedigree.pdf</code>
 */
public class BContainer {

  private final Factor<?>[] mBs;

  /**
   * @param bs value to initialise to
   */
  public BContainer(Factor<?>... bs) {
    mBs = bs;
  }

  /**
   * Set the B value for a given sample within a given family
   * @param familyId unique id within families sample is a parent
   * @param b B value
   */
  public void setB(int familyId, Factor<?> b) {
    mBs[familyId] = b;
  }

  /**
   * Get the B value for a given sample within a given family
   * @param familyId unique id within families sample is a parent
   * @return B value
   */
  public Factor<?> getB(int familyId) {
    return mBs[familyId];
  }

  /**
   * @return number of families in which sample is a parent
   */
  public int size() {
    return mBs.length;
  }
}
