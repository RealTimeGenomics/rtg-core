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
package com.rtg.index.hash;

import java.io.IOException;

import com.rtg.launcher.ISequenceParams;
import com.rtg.mode.Frame;
import com.rtg.reader.SequencesReader;
import com.rtg.util.array.ImmutableIntArray;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.ErrorType;
import com.rtg.util.diagnostic.SlimException;

/**
 */
public abstract class HashLoop {

  protected final int mWindowSize;
  protected final int mStepSize;
  private int mPadding; //TODO this should be final, but will require an inordinate amount of work to force sub classes to comply

  protected HashLoop(int windowSize, int stepSize) {
    mWindowSize = windowSize;
    mStepSize = stepSize;
  }

  /**
   * Implementers should iterate over a set of hashes in the sequences provided and perform a <code>hashCall</code> for each one.
   * It is up to the implementer to determine the set of hashes to generate.
   * @param params parameters for the sequence
   * @param byteBuffer used for reading sequences.
   * @throws IOException if an IO error occurs
   * @return the total number of nucleotides read.
   */
  public abstract long execLoop(ISequenceParams params, byte[] byteBuffer) throws IOException;

  /**
   * This is called once for each complete window.
   * Should be called for each hash in the set generated by <code>execLoop</code>
   * @param hash code generated from the window.
   * @param internalId unique id for each sequence/frame combination.
   * @param stepPosition position in sequence (divided by step size).
   * @throws IOException if an I/O error occurs.
   */
  public abstract void hashCall(final long hash, final int internalId, final int stepPosition) throws IOException;

  /**
   * Hash call for both forward and reverse at once. Applicable to nucleotide only.
   * @param hashForward hash value for forward
   * @param hashReverse hash value for reverse
   * @param stepPosition position in sequence (divided by step size). (for forward)
   * @param internalId unique id for each sequence/frame combination.
   * @throws IOException whenever
   */
  public abstract void hashCallBidirectional(final long hashForward, final long hashReverse, final int stepPosition, final int internalId) throws IOException;


  /**
   * Used for building on paired reads, build loops should override to store which side is being operated on.
   * @param second true if second side
   */
  public void setSide(boolean second) {

  }

  /**
   * Called for each new sequence.
   * @param seqId the internal id of the sequence.
   * @param length of sequence.
   * @throws IOException if an IO error occurs
   */
  public void nextSeq(final int seqId, final int length) throws IOException {
    //default do nothing;
  }

  /**
   * Called for each new sequence/frame pair.
   * @param seq the external id of the sequence.
   * @param frame the frame.
   */
  public abstract void next(final long seq, final Frame frame);

  /**
   * Called after each sequence/frame pair has been completed.
   * @throws IOException if an I/O error occurs.
   */
  public abstract void end() throws IOException;

  /**
   * Called after each sequence has been completed.
   * @throws IOException if an I/O error occurs.
   */
  public abstract void endSequence() throws IOException;

  /**
   * Called at the very end after all sequence/frame pairs have been completed.
   * @throws IOException if error
   */
  public abstract void endAll() throws IOException;

  /**
   * Get the lengths of the reads (indexed on internal ids).
   * @return the lengths of the reads.
   */
  public ImmutableIntArray readLengths() {
    throw new UnsupportedOperationException();
  }


  /**
   * Set the padding on the clip regions for splitting workload
   * @param padding the padding
   */
  public void setThreadPadding(int padding) {
    mPadding = padding;
  }

  /**
   * Get the padding on the clip regions for splitting workload
   * @return the padding
   */
  public int getThreadPadding() {
    return mPadding;
  }

  /**
   * Utility method to create a byte buffer the length of the longest sequence in the supplied reader. This buffer can be used in calls to <code>execLoop</code>.
   *
   * @param reader the reader containing the sequences.
   * @return a byte array of the length of the longest sequence.
   * @throws SlimException if the longest sequence is greater than Integer.MAX_VALUE in length
   */
  public static byte[] makeBuffer(SequencesReader reader) {
    final long maxLength = reader.maxLength();
    if (maxLength >= Integer.MAX_VALUE) {
      Diagnostic.error(ErrorType.SEQUENCE_TOO_LONG, String.valueOf(maxLength));
      throw new SlimException();
    }
    return new byte[(int) maxLength];
  }
}

