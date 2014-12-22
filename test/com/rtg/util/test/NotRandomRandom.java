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
package com.rtg.util.test;

import com.rtg.util.PortableRandom;

/**
 * Not-random random generator for tests
 */
public class NotRandomRandom extends PortableRandom {

  private double mNextDouble;
  private int mNextInt;
  private boolean mNextBoolean;

  /**
   * Constructor
   */
  public NotRandomRandom() {
    super();
    mNextDouble = 0.0;
    mNextInt = 0;
    mNextBoolean = false;
  }

  /**
   */
  @Override
  public double nextDouble() {
    //System.err.println("next double" + mNextDouble );
    final double rand = mNextDouble;
    mNextDouble += 0.1;
    if (mNextDouble >= 1.0) {
      mNextDouble = 0.0;
    }
    //System.err.println("nextDouble: " + rand);
    return rand;
  }

  /**
   */
  @Override
  public int nextInt(int max) {
    if (mNextInt >= max) {
      mNextInt = 0;
    }
    final int rand = mNextInt;
    mNextInt += 1;
    //System.err.println("nextInt: " + rand);
    return rand;
  }

  /**
   */
  @Override
  public boolean nextBoolean() {
    final boolean rand = mNextBoolean;
    mNextBoolean = !mNextBoolean;
    return rand;
  }
}


