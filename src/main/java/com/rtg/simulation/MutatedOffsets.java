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

package com.rtg.simulation;

import java.util.Arrays;

import com.rtg.util.diagnostic.SlimException;

/**
 * A class for managing the arrays for a set of mutation offsets.
 */
public class MutatedOffsets {

  private static final int START_ARRAY_SIZE = 100;

  private static final double ARRAY_GROW_FACTOR = 3.0 / 2.0;

  private int[] mPositions;
  private int[] mOffsets;


  private int mLength = 1;

  private boolean mFrozen = false;

  /**
   * Constructor.
   */
  public MutatedOffsets() {
    this(START_ARRAY_SIZE);
  }

  protected MutatedOffsets(int initialSize) {
    mPositions = new int[initialSize];
    mOffsets = new int[initialSize];
  }

  private void checkFrozen(boolean shouldBeFrozen) {
    if (!shouldBeFrozen && mFrozen) {
      throw new SlimException("Cannot add offsets once frozen.");
    } else if (shouldBeFrozen && !mFrozen) {
      throw new SlimException("Cannot use to get offsets until frozen.");
    }
  }

  /**
   * Get the position in the original reference
   * given a position in the mutated genome.
   * @param position 1-based position in the mutated genome.
   * @return the 1-based position in the original reference.
   */
  public int getOffsetPosition(int position) {
    checkFrozen(true);
    final int index = Arrays.binarySearch(mPositions, position);
    if (index >= 0) {
      return position + mOffsets[index];
    } else {
      return position + mOffsets[Math.abs(index + 1) - 1];
    }
  }

  /**
   * Add an insertion to the offset arrays at the given
   * position of the given size.
   * @param position 1-based position on the original reference (of the anchor base).
   * @param size the size of the insertion.
   */
  public void addInsertion(int position, int size) {
    checkFrozen(false);
    int pos = position - mOffsets[mLength - 1] + 1;
    checkNewPosition(pos);
    for (int i = 0; i <= size; ++i, ++pos) {
      if (i == size / 2) {
        continue;
      }
      addOffset(pos, mOffsets[mLength - 1] - 1);
    }
  }

  /**
   * Add an insertion to the offset arrays at the given
   * position of the given size.
   * @param position 1-based position on the original reference (of the anchor base).
   * @param size the size of the deletion.
   */
  public void addDeletion(int position, int size) {
    checkFrozen(false);
    final int lastOffset = mOffsets[mLength - 1];
    final int pos = position - lastOffset + 1;
    checkNewPosition(pos);
    addOffset(pos, lastOffset + size);
  }

  private void checkNewPosition(int pos) {
    if (pos <= mPositions[mLength - 1]) {
      throw new IllegalArgumentException("Cannot add a position earlier than existing positions.");
    }
  }

  /**
   * Add an offset to the array at the given position.
   * @param position 1-based the position on the mutated genome.
   * @param offset the offset to apply from then on.
   */
  private void addOffset(int position, int offset) {
    if (mPositions.length == mLength) {
      growArrays();
    }
    mPositions[mLength] = position;
    mOffsets[mLength] = offset;
    ++mLength;
  }

  private void growArrays() {
    if (mPositions.length == Integer.MAX_VALUE) {
      throw new SlimException("Too many offsets in VCF to be represented.");
    }
    final int newSize = (int) Math.min(ARRAY_GROW_FACTOR * mPositions.length, Integer.MAX_VALUE);
    mPositions = Arrays.copyOf(mPositions, newSize);
    mOffsets = Arrays.copyOf(mOffsets, newSize);
  }

  /**
   * Freeze the arrays for use in calculating offsets.
   */
  public void freeze() {
    if (!mFrozen) {
      mFrozen = true;
      mPositions = Arrays.copyOf(mPositions, mLength);
      mOffsets = Arrays.copyOf(mOffsets, mLength);
    }
  }
}
