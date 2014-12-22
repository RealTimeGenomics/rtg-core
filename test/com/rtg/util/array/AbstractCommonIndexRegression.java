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

package com.rtg.util.array;

import java.util.Random;

import com.rtg.util.diagnostic.Diagnostic;

import junit.framework.TestCase;

/**
 * Base class for testing all of the index implementations that are capable of
 * holding more than the maximum integer worth of data that they actually work.
 */
public abstract class AbstractCommonIndexRegression extends TestCase {

  private static final long NUM_ELEMENTS = 2L * Integer.MAX_VALUE + 9000L;

  protected abstract long getRange();
  protected abstract CommonIndex createIndex(long elements);

  protected long getNumElements() {
    return NUM_ELEMENTS;
  }

  /**
   * Test the common index implementation
   */
  public void testIndex() {
    Diagnostic.setLogStream();
    doTest(getRange(), getNumElements());
  }

  private void doTest(long range, long elements) {
    final RandomLongGenerator value = new RandomLongGenerator(range);

    final CommonIndex index = createIndex(elements);
    assertEquals(elements, index.length());

    for (long l = 0; l < elements; l++) {
      index.set(l, value.nextValue());
    }

    value.reset();

    for (long l = 0; l < elements; l++) {
      assertEquals(value.nextValue(), index.get(l));
    }
  }

  /**
   * Utility class to create a repeatable stream of
   * positive random longs.
   */
  public static final class RandomLongGenerator {

    private final long mRange;
    private final long mSeed;
    private Random mRand;

    /**
     * Constructor for generator with seed from current time.
     * @param range the range of longs to generate from 0 inclusive to range exclusive.
     */
    public RandomLongGenerator(long range) {
      this(range, System.nanoTime());
      System.err.println("Seed: " + mSeed);
    }

    /**
     * Constructor for generator with seed provided.
     * @param range the range of longs to generate from 0 inclusive to range exclusive.
     * @param seed the seed for the random number generator.
     */
    public RandomLongGenerator(long range, long seed) {
      assert 0 < range && range <= Long.MAX_VALUE;
      mRange = range;
      mSeed = seed;
      reset();
    }

    /**
     * Method to get the next long in the sequence.
     * @return the next long in the sequence.
     */
    public long nextValue() {
      return (long) (mRand.nextDouble() * mRange);
    }

    /**
     * Method to reset the generator to produce the same sequence again
     * from the beginning.
     */
    public void reset() {
      mRand = new Random(mSeed);
    }
  }
}
