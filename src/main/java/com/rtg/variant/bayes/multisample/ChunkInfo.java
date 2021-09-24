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
