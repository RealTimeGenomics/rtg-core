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
package com.rtg.blacklist;

import com.rtg.util.MathUtils;
import com.rtg.util.array.atomic.AtomicIndex;
import com.rtg.util.array.atomic.AtomicIntChunks;
import com.rtg.util.array.atomic.AtomicLongChunks;


/**
 * Lock free thread safe (for increment method only!) hash counter
 */
public class HashCounter {

  private static final int BOUNCE_LIMIT = 126;
  private static final int BOUNCE_BITS = MathUtils.ceilPowerOf2Bits(BOUNCE_LIMIT);
  private static final int[] BOUNCE = new int[BOUNCE_LIMIT];
  static {
    for (int i = 0; i < BOUNCE_LIMIT; i++) {
      BOUNCE[i] = i * (i + 1) / 2;
    }
  }

  //upper_hash_bits|num_bounces|hit_count
  private final AtomicIndex mHashes;
  private final long mPositionMask; //turn shuffled key into position in hash array
  private final int mMaskBits; //bits for above mask
  private final long mKeyMask; //masks out the stored (upper) bits of the key
  private final long mBounceMask; //masks out the number of bounces
  private final long mLength; //length of mHashes
  private final BinaryMatrix mMatrix; //shuffles (randomises) key in reversible fashion
  private final BinaryMatrix mReverseMatrix; //reverses the shuffle
  private final int mCountBits; //number of bits for the counter
  private final long mCountMask; //mask for the counter

  //iterator
  private long mCurrentIndex = -1;
  private long mCurrentKey;
  private long mCurrentCount;


  /**
   * Creates a hash counter
   * @param length minimum number of hash entries
   * @param keyBits number of bits required to represent key
   * @param maxRepeat maximum number required to be stored in count field
   */
  public HashCounter(long length, int keyBits, int maxRepeat) {
    mHashes = minMemFor(length, keyBits, maxRepeat);
    mLength = mHashes.length();
    mPositionMask = mLength - 1;
    mMaskBits = MathUtils.ceilPowerOf2Bits(mPositionMask);
    mKeyMask = (1L << (keyBits - mMaskBits)) - 1;
    mBounceMask = (1L << BOUNCE_BITS) - 1;
    final int storageBits = mHashes.fieldBits();
    mCountBits = storageBits - (keyBits - mMaskBits) - BOUNCE_BITS;
    mCountMask = (1L << mCountBits) - 1;
    mMatrix = BinaryMatrix.createReversibleMatrix(keyBits);
    mReverseMatrix = mMatrix.invert();
  }

  /**
   * Creates the smallest (in memory use) array for storing the hashes that meets the following criteria
   * @param length minimum number of hash entries
   * @param keyBits number of bits required to represent key
   * @param maxRepeat maximum number required to be stored in count field
   * @return the array
   */
  static AtomicIndex minMemFor(long length, int keyBits, int maxRepeat) {
    final int minPositionBits = MathUtils.ceilPowerOf2Bits(length);
    final int minCountBits = MathUtils.ceilPowerOf2Bits(maxRepeat);

    //mimimum space for intarray
    final int positionBitsInt = Math.max(keyBits + minCountBits + BOUNCE_BITS - 32, minPositionBits);
    //minimum space for longarray
    final int positionBitsLong = Math.max(keyBits + minCountBits + BOUNCE_BITS - 64, minPositionBits);
    //find the options which uses least memory
    if (positionBitsLong <= positionBitsInt - 1) { //using 1 extra bit for length in ints will be same amount of memory for longs
      final long actualLength = 1L << positionBitsLong;
      return new AtomicLongChunks(actualLength);
    } else {
      final long actualLength = 1L << positionBitsInt;
      return new AtomicIntChunks(actualLength);
    }
  }

  /**
   * increments the count for given key (this method is thread safe).
   * @param key the key
   * @throws TooManyCollisionsException if the key bounces too many times whilst looking for a location
   */
  public void increment(long key) throws TooManyCollisionsException {
    final long hash = mMatrix.times(key);
    final long hashUpperBits = (hash >>> mMaskBits) & mKeyMask;
    final long originalPos = hash & mPositionMask;
    long bounces = 0;
    long pos = -1;
    long getVal = 0; //logically the value currently stored at pos
    while (bounces < BOUNCE_LIMIT) {
      pos = (originalPos + BOUNCE[(int) bounces]) & mPositionMask;
      getVal = mHashes.get(pos);
      //is position empty?
      final long getBounce = (getVal >>> mCountBits) & mBounceMask;
      if (getBounce == 0) { //bounce field == 0 implies unclaimed position
        //claim spot
        final long setVal = hashUpperBits << (mCountBits + BOUNCE_BITS)
                      | ((bounces + 1) << mCountBits);
        if (mHashes.compareAndSet(pos, 0, setVal)) {
          getVal = setVal;
          break;
        }
        //something else claimed our spot
        //reloop
      } else {
        //is it this hash?
        if (bounces == getBounce - 1) { //same bounce value (since we store bounces + 1)
          final long getUpperBits = (getVal >>> (mCountBits + BOUNCE_BITS)) & mKeyMask;
          if (getUpperBits == hashUpperBits) { //same upper bits
            //same hash
            break;
          }
        }
        //different hash
        bounces++;
        //reloop
      }
    }
    if (bounces >= BOUNCE_LIMIT) {
      throw new TooManyCollisionsException("Exceed maximum number of bounces in hash table");
    }
    //pos should now have position of hash and getVal has current value at pos
    do {
      final long count = getVal & mCountMask;
      if (count == mCountMask) {
        break; //maximum value
      }
      if (mHashes.compareAndSet(pos, getVal, getVal + 1)) {
        //incremented
        break;
      }
      //failed;
      getVal = mHashes.get(pos);
    } while (true);
  }

  private long posToKey(long pos, long upperBits, long bounce) {
    return (upperBits << mMaskBits) | ((pos - BOUNCE[(int) bounce]) & mPositionMask);
  }

  /**
   * Sets {@link HashCounter#getKey()} or {@link HashCounter#getCount()} to the next value
   * All calls to increment should be finished before any calls to this.
   * @return true if there is another value
   */
  public boolean next() {
    while (mCurrentIndex + 1 < mLength) {
      mCurrentIndex++;
      final long currentVal = mHashes.get(mCurrentIndex);
      if (currentVal != 0) {
        mCurrentCount = currentVal & mCountMask;
        long tmp = currentVal >>> mCountBits;
        final int bounce = (int) (tmp & mBounceMask) - 1;
        tmp = tmp >>> BOUNCE_BITS;
        mCurrentKey = posToKey(mCurrentIndex, tmp & mKeyMask, bounce);
        return true;
      }
    }
    return false;
  }

  /**
   * gets the current key
   * @return the key
   */
  public long getKey() {
    return mReverseMatrix.times(mCurrentKey);
  }

  /**
   * gets the current count
   * @return the count
   */
  public long getCount() {
    return mCurrentCount;
  }

  /**
   * Indicates that the hash map failed too many times for a single lookup
   */
  public static class TooManyCollisionsException extends Exception {
    /**
     * @param message  see {@link Exception}
     */
    public TooManyCollisionsException(String message) {
      super(message);
    }
  }

}
