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

import java.util.Random;

/**
 * Utility class to create a repeatable stream of
 * positive random bytes, maximum range 128 (0 - 127).
 */
public final class RandomByteGenerator {

  private final int mRange;
  private final long mSeed;
  private Random mRand;

  /**
   * Constructor for generator with seed from current time.
   * @param range the range of bytes to generate from 0 inclusive to range exclusive.
   */
  public RandomByteGenerator(int range) {
    this(range, System.nanoTime());
    System.err.println("Seed: " + mSeed);
  }

  /**
   * Constructor for generator with seed provided.
   * @param range the range of bytes to generate from 0 inclusive to range exclusive.
   * @param seed the seed for the random number generator.
   */
  public RandomByteGenerator(int range, long seed) {
    assert 0 < range && range <= 128;
    mRange = range;
    mSeed = seed;
    reset();
  }

  /**
   * Method to get the next byte in the sequence.
   * @return the next byte in the sequence.
   */
  public byte nextValue() {
    return (byte) mRand.nextInt(mRange);
  }

  /**
   * Method to reset the generator to produce the same sequence again
   * from the beginning.
   */
  public void reset() {
    mRand = new Random(mSeed);
  }
}
