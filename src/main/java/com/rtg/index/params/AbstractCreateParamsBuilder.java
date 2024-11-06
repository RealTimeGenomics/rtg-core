/*
 * Copyright (c) 2018. Real Time Genomics Limited.
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
package com.rtg.index.params;

import com.reeltwo.jumble.annotations.TestClass;

/**
 * Abstract base builder for create params
 */
@TestClass("com.rtg.index.params.CreateParamsTest")
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

  protected abstract B self();
}
