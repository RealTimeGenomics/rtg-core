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
package com.rtg.util;

/**
 * Portable random number generator for cross-compatibility with Java and .NET
 * Based on the Knuth algorithm, using the java values
 * @see <a href="http://www.particle.kth.se/~lindsey/JavaCourse/Book/Part1/Tech/Chapter04/javaRandNums.html">Java Rand Nums</a>
 */
public class PortableRandom {

  private static final long MULTIPLIER = 0x5DEECE66DL;
  private static final long ADDEND = 0xBL;
  private static final long MASK = (1L << 48) - 1;
  private static long sSeedUniquifier = 8682522807148012L;
  private static synchronized long nextSeedUniquifier() {
    return ++sSeedUniquifier;
  }

  private double mNextNextGaussian;
  private boolean mHaveNextNextGaussian = false;
  private long mSeed;
  private long mInitialSeed;

  /**
   * Create a new Portable random number generator seeded by time
   */
  public PortableRandom() {
    this(nextSeedUniquifier() + System.nanoTime());
  }

  /**
   * Create a new Portable random number generator with a specified seed
   * @param seed the seed
   */
  public PortableRandom(long seed) {
    setSeed(seed);
  }

  protected int next(int bits) {
    mSeed = (mSeed * MULTIPLIER + ADDEND) & MASK;
    return (int) (mSeed >>> (48 - bits));
  }

  /**
   * Reset the state of this PortableRandom using the given seed.
   * @param seed the seed
   */
  public void setSeed(long seed) {
    mInitialSeed = seed;
    mSeed = (mInitialSeed ^ MULTIPLIER) & MASK;
    mHaveNextNextGaussian = false;
  }

  /**
   * @return the next random boolean
   */
  public boolean nextBoolean() {
    return next(1) != 0;
  }

  /**
   * @return the next random double between 0.0 (inclusive) and 1.0 (exclusive)
   */
  public double nextDouble() {
    return (((long) next(26) << 27) + next(27)) / (double) (1L << 53);
  }

  /**
   * @return the next random float
   */
  //public float nextFloat() {
  //  return next(24) / ((float) (1 << 24));
  //}

  /**
   * @return the next random integer
   */
  public int nextInt() {
    return next(32);
  }

  /**
   * @return the next random long
   */
  public long nextLong() {
    return ((long) next(32) << 32) + next(32);
  }


  /**
   * Return the next Gaussian
   * See Knuth, ACP, Section 3.4.1 Algorithm C.
   * @return the next Gaussian double
   */
  public double nextGaussian() {
    if (mHaveNextNextGaussian) {
      mHaveNextNextGaussian = false;
      return mNextNextGaussian;
    } else {
      double v1, v2, s;
      do {
        v1 = 2 * nextDouble() - 1; // between -1 and 1
        v2 = 2 * nextDouble() - 1; // between -1 and 1
        s = v1 * v1 + v2 * v2;
      } while (s >= 1 || s == 0);
      final double multiplier = StrictMath.sqrt(-2 * StrictMath.log(s) / s);
      mNextNextGaussian = v2 * multiplier;
      mHaveNextNextGaussian = true;
      return v1 * multiplier;
    }
  }

  private static final int INTEGER_SIZE = 32;
  private static final int BYTE_SIZE = 8;

  /**
   * Fill a specified byte array with random bytes
   * @param bytes the array to fill with random bytes
   */
  public void nextBytes(byte[] bytes) {
    final int length = bytes.length;
    for (int i = 0; i < length;) {
      for (int rnd = nextInt(), n = Math.min(length - i, INTEGER_SIZE / BYTE_SIZE); n-- > 0; rnd >>= BYTE_SIZE) {
        bytes[i++] = (byte) rnd;
      }
    }
  }

  /**
   * Return a random integer between 0 (inclusive) and n (exclusive)
   * @param exclMaxBound exclusive maximum bound for returned random integer
   * @return random integer
   */
  public int nextInt(int exclMaxBound) {
    if (exclMaxBound <= 0) {
      throw new IllegalArgumentException("n must be positive");
    }

    if ((exclMaxBound & -exclMaxBound) == exclMaxBound) {
      return (int) ((exclMaxBound * (long) next(31)) >> 31);
    }
    int bits, val;
    do {
      bits = next(31);
      val = bits % exclMaxBound;
    } while (bits - val + (exclMaxBound - 1) < 0);
    return val;
  }


  /**
   * Method to get the seed value for logging purposes.
   * @return the current initial seed value
   */
  public long getSeed() {
    return mInitialSeed;
  }
}
