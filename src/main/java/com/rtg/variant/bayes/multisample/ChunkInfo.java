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

package com.rtg.variant.bayes.multisample;

/**
 * Shared information about chunks etc.
 */
public class ChunkInfo {

  private final int mChunkSize;
  private final int mZeroBasedStart;
  private final int mZeroBasedEnd;
  private final int mNumberChunks;
  private final int mBufferSize;
  private final int mLength;
  private final String mRefName;

  /**
   * @param seqLength reference sequence length
   * @param refName reference sequence name
   * @param chunkSize size of the chunk
   * @param restrictionStart the start position (zero based, inclusive) of restricted calling. -1 means no position filtering.
   * @param restrictionEnd the end position (zero based, exclusive) of restricted calling. -1 means no position filtering.
   * @param execThreads number of execution threads
   * @param maxReadLength maximum read length
   */
  public ChunkInfo(int seqLength, String refName, int chunkSize, int restrictionStart, int restrictionEnd, int execThreads, int maxReadLength) {
    mChunkSize = chunkSize;
    // Normally we start at nt 0, but if the user has specified a restriction
    // use that instead.
    mZeroBasedStart = restrictionStart < 0 ? 0 : restrictionStart;
    mZeroBasedEnd = restrictionEnd < 0 ? seqLength : restrictionEnd;
    assert mZeroBasedEnd >= mZeroBasedStart;
    mLength = mZeroBasedEnd - mZeroBasedStart;
    mNumberChunks = ((mZeroBasedEnd - mZeroBasedStart) + mChunkSize - 1) / mChunkSize;
    mBufferSize = mChunkSize * (execThreads + 2) * 2 + maxReadLength;
    mRefName = refName;
  }

  /**
   * @return chunk size
   */
  public int chunkSize() {
    return mChunkSize;
  }

  /**
   * @return start position (0-based, inclusive)
   */
  public int start() {
    return mZeroBasedStart;
  }

  /**
   * @return end position (0-based, exclusive)
   */
  public int end() {
    return mZeroBasedEnd;
  }

  /**
   * @return length of region being chunked
   */
  public int length() {
    return mLength;
  }

  /**
   * @param pos position of last completed base
   * @return the percentage of completion
   */
  public int percent(final int pos) {
    return (int) ((100L * (pos - start())) / length());
  }

  /**
   * @return the number of chunks
   */
  public int numberChunks() {
    return mNumberChunks;
  }

  int bufferSize() {
    return mBufferSize;
  }

  /**
   * @return reference sequence name
   */
  public String refName() {
    return mRefName;
  }
}
