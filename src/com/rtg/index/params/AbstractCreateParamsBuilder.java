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

import com.rtg.index.Index;

/**
 * @author Dave Ware
 */
public abstract class AbstractCreateParamsBuilder<B extends AbstractCreateParamsBuilder<B>> {
  protected long mSize = 20;
  protected int mHashBits = 20;
  protected int mWindowBits = 20;
  protected int mValueBits = 20;
  protected boolean mCompressHashes = true;
  protected boolean mCreateBitVector = true;
  protected boolean mSpaceEfficientButUnsafe = false;
  protected boolean mIdeal = false;
  protected boolean mOnlyKeepRepeatHashes = false;
  protected Index mHashFilterIndex = null;

  /**
   * @param size upper bound of number of hash windows expected in index.
   * @return self() builder for chaining purposes
   */
  public B size(long size) {
    mSize = size;
    return self();
  }

  /**
   *
   * @param hashBits number of bits recorded in each hash window.
   * @return self() builder for chaining purposes
   */
  public B hashBits(int hashBits) {
    mHashBits = hashBits;
    return self();
  }

  /**
   * @param windowBits number of bits needed to uniquely represent a window (the hash may lose information).
   * @return self() builder for chaining purposes
   */
  public B windowBits(int windowBits) {
    mWindowBits = windowBits;
    return self();
  }

  /**
   * @param valueBits the number of bits needed to represent the values associated with each hash.
   * @return self() builder for chaining purposes
   */
  public B valueBits(int valueBits) {
    mValueBits = valueBits;
    return self();
  }

  /**
   * @param compressHashes if we are compressing the hash array (requires two passes of adds). <code>noOverflow</code> must be on for this to be allowed
   * @return self() builder for chaining purposes
   */
  public B compressHashes(boolean compressHashes) {
    mCompressHashes = compressHashes;
    return self();
  }

  /**
   * @param createBitVector if we are creating the bit vector
   * @return self() builder for chaining purposes
   */
  public B createBitVector(boolean createBitVector) {
    mCreateBitVector = createBitVector;
    return self();
  }

  /**
   * Implies that the hashes for the index are packed as small as possible, with the tradeoff
   * that potentially an index storage implementation that <b>is not</b> safe from word tearing might
   * be selected. This means it will not be safe to do a multi-threaded sort of the hashes.
   * @param fanatical if we want to use the <code>spaceEfficientButUnsafe</code> <code>ArrayType</code> selection mode
   * @return self() builder for chaining purposes
   */
  public B spaceEfficientButUnsafe(boolean fanatical) {
    mSpaceEfficientButUnsafe = fanatical;
    return self();
  }

  /**
   * Sets whether to minimise the memory used.
   * Only use if not intending to search the index.
   * @param ideal if true then minimise total of hash and initial pointer memory
   * @return self() builder for chaining purposes
   */
  public B ideal(boolean ideal) {
    mIdeal = ideal;
    return self();
  }

  /**
   * @param val true to discard all hashes except those that exceed repeat threshold
   * @return self() builder for chaining purposes
   */
  public B onlyKeepRepeatHashes(boolean val) {
    mOnlyKeepRepeatHashes = val;
    return self();
  }

  /**
   * @param index index that contains hashes which should be discarded
   * @return self() builder for chaining purposes
   */
  public B hashFilterIndex(Index index) {
    mHashFilterIndex = index;
    return self();
  }

  protected abstract B self();
}
