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
package com.rtg.index.params;


import java.io.IOException;

import com.rtg.index.HashBitHandle;
import com.rtg.launcher.BuildParams;
import com.rtg.util.MathUtils;
import com.rtg.util.StringUtils;
import com.rtg.util.Utils;
import com.rtg.util.array.ArrayHandle;
import com.rtg.util.array.ArrayType;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.Integrity;

/**
 * Holds the minimal set of parameters that are sufficient to enable the initial creation
 * of an <code>IndexImplementation</code> and calculation of the resulting memory usage.
 *
 * Also provides estimates of the amount of memory to be used.
 *
 * <p>
 * WARNING: this class is very intimately tied with the <code>IndexImplementation</code>
 * class. This is required by the need to be able to compute the memory size of
 * <code>IndexImplementation</code> without actually creating it.
 */
public class CreateParams implements Integrity {

  private static final int LONG_BITS = 64;

  /** The multiplier for the length of the bit Vector. */
  public static final long SPEED_UP_FACTOR = 15;

  //Total number of hash bits before any compression
  private final int mHashBits;

  private final int mWindowBits;

  private final int mHashCompressedBits;

  private final int mInitialPointerBits;

  private final long mSize;

  private final int mValueBits;

  private final boolean mCompressHashes;

  private final boolean mCreateBitVector;

  private final boolean mSpaceEfficientButUnsafe;

  private final boolean mOnlyKeepRepeatHashes;

  /**
   * @param size upper bound of number of hash windows expected in index.
   * @param hashBits number of bits recorded in each hash window.
   * @param windowBits number of bits needed to uniquely represent a window (the hash may lose information).
   * @param valueBits the number of bits needed to represent the values associated with each hash.
   * @param compressHashes if we are compressing the hash array (requires two passes of adds). <code>noOverflow</code> must be on for this to be allowed
   * @param createBitVector if we are creating the bit vector;
   * @param spaceEfficientButUnsafe if we want to use the <code>spaceEfficientButUnsafe</code> <code>ArrayType</code> selection mode
   * @param ideal if true then minimise total of hash and initial pointer memory
   */
  public CreateParams(final long size, final int hashBits, final int windowBits, final int valueBits, final boolean compressHashes, boolean createBitVector, boolean spaceEfficientButUnsafe, boolean ideal) {
    //System.err.println("size=" + size + " hashBits=" + hashBits + " windowBits=" + windowBits + " valueBits=" + valueBits);
    mSize = size;
    mHashBits = hashBits;
    mWindowBits = windowBits;
    mValueBits = valueBits;
    mCompressHashes = compressHashes;
    mCreateBitVector = createBitVector;
    mSpaceEfficientButUnsafe = spaceEfficientButUnsafe;
    if (ideal) {
      mInitialPointerBits = computeInitialPointerBitsIdeal(mSize, mHashBits, mSpaceEfficientButUnsafe);
    } else {
      mInitialPointerBits = computeInitialPointerBits(mHashBits, mSize);
    }
    mHashCompressedBits = computeHashCompressedBits(mHashBits, mInitialPointerBits, compressHashes);
    mOnlyKeepRepeatHashes = false;
    //System.err.println(this);
    assert localIntegrity();
  }

  /**
   * @param builder builder containing parameter values
   */
  public CreateParams(final AbstractCreateParamsBuilder<?> builder) {
    mSize = builder.mSize;
    mHashBits = builder.mHashBits;
    mWindowBits = builder.mWindowBits;
    mValueBits = builder.mValueBits;
    mCompressHashes = builder.mCompressHashes;
    mCreateBitVector = builder.mCreateBitVector;
    mSpaceEfficientButUnsafe = builder.mSpaceEfficientButUnsafe;
    if (builder.mIdeal) {
      mInitialPointerBits = computeInitialPointerBitsIdeal(mSize, mHashBits, mSpaceEfficientButUnsafe);
    } else {
      mInitialPointerBits = computeInitialPointerBits(mHashBits, mSize);
    }
    mHashCompressedBits = computeHashCompressedBits(mHashBits, mInitialPointerBits, builder.mCompressHashes);
    mOnlyKeepRepeatHashes = builder.mOnlyKeepRepeatHashes;
    assert localIntegrity();
  }

  /**
   * Compute the initial pointer size that gives the minimum total memory for hash table and initial pointer table.
   * Do not use this if you intend to search an index.
   * @param size total number of entries
   * @param hashBits number of hash bits before adjusting for compression
   * @param fanatical if true use bit packing for hash array
   * @return the number of bits in the initial pointer
   */
  static int computeInitialPointerBitsIdeal(long size, int hashBits, boolean fanatical) {
    //System.err.println("size:" + size + " hashBits:" + hashBits);
    final int def = computeInitialPointerBits(hashBits, size);
    //System.err.println("def:" + def);
    final ArrayType ip = ArrayType.bestForLength(size);
    //System.err.println("ip:" + ip);
    if (def <= 9) {
      return def;
    }
    long bestBytes = Long.MAX_VALUE;
    int best = def;
    for (int l = Math.min(def, hashBits); l >= 10 && hashBits - l <= 64; --l) {
      final long ipBytes = ip.bytes(1L << l);
      final long hashBytes;
      if (fanatical) {
        hashBytes = ArrayType.bestForBitsSpaceEfficientButNotSafeFromWordTearing(hashBits - l).bytes(size);
      } else {
        hashBytes = ArrayType.bestForBitsAndSafeFromWordTearing(hashBits - l).bytes(size);
      }
      if (ipBytes < 0 || hashBytes < 0) {
        //deal with extreme cases where even longs overflow - makes this code more robust
        continue;
      }
      final long totalBytes = ipBytes + hashBytes;
      //System.err.println("l:" + l + " ipBytes:" + ipBytes + " hashBytes:" + hashBytes + " totalBytes:" + totalBytes);
      if (totalBytes < bestBytes) {
        bestBytes = totalBytes;
        best = l;
      }
    }
    return best;
  }

  private static int computeInitialPointerBits(int hashBits, long size) {
    //Try for a size which rounds length up to a power of 2 and then decreases by 2 (giving at worst a density of 4)
    final int pBits;
    final int nbits = MathUtils.ceilPowerOf2Bits(size) - 1;
    if (nbits <= 2) { //avoid tricky edge cases with signs etc.
      pBits = 2;
    } else {
      pBits = nbits;
    }
    //make sure big enough to fit the compressed hashes into a long
    final int pointerBits;
    if (hashBits - pBits > Long.SIZE) {
      pointerBits = hashBits - Long.SIZE;
    } else {
      pointerBits = pBits;
    }
    return pointerBits;
  }

  private static int computeHashCompressedBits(int hashBits, int initialPointerBits, boolean compressHashes) {
    final int reduction = compressHashes ? initialPointerBits : 0;
    final int bits0 = hashBits - reduction;
    return Math.max(0, bits0);
  }

  /**
   * @return an upper bound on the number of hash windows to be generated
   * by the subjects.
   */
  public long size() {
    return mSize;
  }

  /**
   * Compute the number of bits required to store a hash code.
   * @return the number of bits required to store a hash code.
   */
  public int hashBits() {
    return mHashBits;
  }

  /**
   * @return the number of bits used in hash after allowing for compression.
   */
  public int hashCompressedBits() {
    return mHashCompressedBits;
  }

  /**
   * Get the number of bits required to uniquely represent a window.
   * @return the number of bits required to uniquely represent a window.
   */
  public int windowBits() {
    return mWindowBits;
  }

  /**
   * @return true if we are compressing the hash array (requires two passes of adds)
   */
  public boolean compressHashes() {
    return mCompressHashes;
  }

  private ArrayType bestForBits(final int bits) {
    if (mSpaceEfficientButUnsafe) {
      return ArrayType.bestForBitsSpaceEfficientButNotSafeFromWordTearing(bits);
    } else {
      return ArrayType.bestForBitsAndSafeFromWordTearing(bits);
    }
  }

  /**
   * Construct a handle which can be used to construct a hash array for <code>IndexImplementation</code>.
   * @return a handle which can be used to construct a hash array for <code>IndexImplementation</code>.
   */
  public ArrayHandle hash() {
    final int bits = mHashCompressedBits;
    final ArrayType type = bestForBits(bits >= 0 ? bits : 0);
    //System.err.println("hash twoPass()=" + twoPass() + " initialPointerBits()=" + initialPointerBits() + " mHashBits=" + mHashBits + " bits=" + bits + " compact()=" + compact() + " type=" + type);
    return new ArrayHandle(type, mSize);
  }

  /**
   * Get the number of bits required to uniquely represent the value associated with
   * a hash code..
   * @return the number of bits required to uniquely represent a value.
   */
  public int valueBits() {
 //   new RuntimeException("I just want a stack trace damnit: " + mValueBits).printStackTrace();
    return mValueBits;
  }

  /**
   * Construct a handle which can be used to construct a a value array for <code>IndexImplementation</code>.
   * @return a handle which can be used to construct a value array for <code>IndexImplementation</code>.
   */
  public ArrayHandle value() {
    final ArrayType type = bestForBits(valueBits());
    //System.err.println("value valueBits()=" + valueBits() + " type=" + type);
    return new ArrayHandle(type, mSize);
  }

  /**
   * Compute the number of bits used to compute the length of the initial position array.
   * @return the number of bits used to compute the length of the initial position array.
   */
  public int initialPointerBits() {
    return mInitialPointerBits;
  }

  /**
   * Construct a handle which can be used to construct an initial position array for <code>IndexImplementation</code>.
   * @return a handle which can be used to construct an initial position array for <code>IndexImplementation</code>.
   */
  public ArrayHandle initialPosition() {
    //extra entries at end to simplify code
    //first +1 allows ip[i]...ip[i+1] to be used as start to exclusive end when pulling out regions
    //second +1 allows for the entire array to be shifted right by 1 during first pass of two-pass compressed memory
    final long initSize =  (1L << initialPointerBits()) + 2;
    final ArrayType type = ArrayType.bestForLength(mSize + 1);
    return new ArrayHandle(type, initSize);
  }

  /**
   * Construct a handle which can be used to construct a bit vector for <code>IndexImplementation</code>.
   * @return a handle which can be used to construct a bit vector for <code>IndexImplementation</code>.
   */
  public HashBitHandle bitVector() {
    if (mCreateBitVector) {
      //The bit vector must be large enough to contain 1 entry for each of mSize entries
      //Beyond that it is scaled so that it is speed up factor times mSize except that there is no point in
      //exceeding mHashBits as then there will be 1 bit per hash.
      final int ms = MathUtils.ceilPowerOf2Bits(mSize);
      final int sz = MathUtils.ceilPowerOf2Bits(mSize * SPEED_UP_FACTOR);
      final int bv = Math.max(ms, Math.min(sz, mHashBits));
      //System.err.println("ms=" + ms + " sz=" + sz + " bv=" + bv + " mSize=" + mSize + " SPEED_UP_FACTOR=" + SPEED_UP_FACTOR);
      assert bv >= ms;
      return new HashBitHandle(Math.min(Long.SIZE, mHashBits), bv);
    } else {
      return null;
    }
  }

  @Override
  public int hashCode() {
    return Utils.pairHash(mHashBits, hashLong(mSize));
  }

  private int hashLong(final long l) {
    return (int) l ^ (int) (l >> 32);
  }

  @Override
  public boolean equals(final Object obj) {
    if (this == obj) {
      return true;
    }
    if (obj == null) {
      return false;
    }
    final CreateParams that = (CreateParams) obj;
    return this.mHashBits == that.mHashBits
        && this.mSize == that.mSize;
  }


  @Override
  public boolean integrity() {
    localIntegrity();
    return true;
  }

  /**
   * Class invariant. Returns true so can be used in assert statements.
   * This is here so it can be called from the constructor without
   * calling the overridden integrity method in the child classes.
   */
  private boolean localIntegrity() {
    Exam.assertTrue("" + size(), size() >= 0);
    Exam.assertTrue(hashBits() >= 1); // && hashBits() <= LONG_BITS);
    //Exam.assertTrue(hashCompressedBits() <= LONG_BITS);
    if (windowBits() <= LONG_BITS) {
      Exam.assertTrue(hashBits() == windowBits());
    } else {
      Exam.assertTrue(hashBits() == LONG_BITS || hashBits() == windowBits());
    }
    Exam.assertTrue(initialPointerBits() >= 2);
//    Exam.assertTrue(valueBits() >= 0 && valueBits() <= LONG_BITS);
    return true;
  }

  @Override
  public String toString() {
    return " size=" + StringUtils.commas(mSize)
        + " hash bits=" + mHashBits
        + " initial pointer bits=" + initialPointerBits()
        + " value bits=" + valueBits();
  }

  @Override
  public boolean globalIntegrity() {
    integrity();
    return true;
  }


  /**
   * Calculate the number of hash bits needed.
   * @param typeBits bits needed for each residue.
   * @param windowSize number of residues in a hash window.
   * @return the number of required hash bits.
   */
  public static int calculateHashBits(final int typeBits, final int windowSize) {
    final int winBits = windowSize * typeBits;
    return winBits > 64 ? 64 : winBits;
  }

  /**
   * builder for create params with defaults calculated from build params
   * @param bp the build params
   * @return builder with some of values preset
   * @throws IOException if an IO error occurs whilst calculating size
   */
  public static CreateParamsBuilder fromBuildParams(BuildParams bp) throws IOException {
    final int hashBits = CreateParams.calculateHashBits(bp.sequences().mode().codeType().bits(), bp.windowSize());
    final CreateParamsBuilder ret = new CreateParamsBuilder();
    ret.windowBits(hashBits)
      .size(bp.calculateSize())
      .hashBits(hashBits)
      .spaceEfficientButUnsafe(false)
      .ideal(false);
    return ret;
  }


  /**
   * Builder class for <code>CreateParams</code>
   */
  public static class CreateParamsBuilder extends AbstractCreateParamsBuilder<CreateParamsBuilder> {

    /**
     * @return a CreateParams object as described by the builder
     */
    public CreateParams create() {
      return new CreateParams(mSize, mHashBits, mWindowBits, mValueBits, mCompressHashes, mCreateBitVector, mSpaceEfficientButUnsafe, mIdeal);
    }

    @Override
    protected CreateParamsBuilder self() {
      return this;
    }
  }

}
