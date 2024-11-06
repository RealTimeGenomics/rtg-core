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
package com.rtg.launcher;

import java.io.File;
import java.io.IOException;

import com.rtg.mode.SequenceMode;
import com.rtg.reader.SequencesReader;
import com.rtg.reference.Sex;
import com.rtg.util.intervals.LongRange;
import com.rtg.util.Utils;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;


/**
 */
public class MockSequenceParams extends IntegralAbstract implements ISequenceParams {

  private final ReaderParams mReaderParams;
  private final SequenceMode mMode;
  private final HashingRegion mRegion;
  private final LongRange mReaderRestriction;

  /**
   * Constructor that takes reader directly.
   * Only used in testing.
   * @param sr Reader to get sequences from
   * @param mode sequence mode.
   */
  public MockSequenceParams(final SequencesReader sr, final SequenceMode mode) {
    mReaderParams = new MockReaderParams(sr);
    mRegion = new HashingRegion(0, sr.numberSequences());
    mReaderRestriction = LongRange.NONE;
    mMode = mode;
    //integrity();
  }

  /**
   * @param readerParams specifies the directory and provides a reader when needed.
   * @param mode sequence mode.
   */
  public MockSequenceParams(final ReaderParams readerParams, SequenceMode mode) {
    this(readerParams, mode, 0, readerParams.reader().numberSequences());
  }

  /**
   * @param readerParams specifies the directory and provides a reader when needed.
   * @param mode sequence mode.
   * @param start the first sequence to be processed.
   * @param end one past the last sequence to be processed.
   */
  public MockSequenceParams(final ReaderParams readerParams, final SequenceMode mode, final long start, final long end) {
    mReaderParams = readerParams;
    mMode = mode;
    mRegion = new HashingRegion(start, end);
    mReaderRestriction = LongRange.NONE;
    integrity();
  }

  /**
   * Create a new version of this object whose start and end positions lie within the current range.
   * @param region the region to generate params for
   * @return the new <code>SequenceParams</code>
   */
  @Override
  public ISequenceParams subSequence(HashingRegion region) {
    //System.err.println("MockSequenceParams duplicate bingo");
    if (region.getStart() < mRegion.getStart() || region.getEnd() > mRegion.getEnd()) {
      throw new RuntimeException("start=" + region.getStart() + " end=" + region.getEnd() + " this=" + this);
    }
    //Warning: mReaderParams is unsafe in multi-core - not clear if some form of duplication may be needed for testing
    return new MockSequenceParams(mReaderParams, mode(), region.getStart(), region.getEnd());
  }

  /**
   * Get the mode.
   * @return the mode.
   */
  @Override
  public SequenceMode mode() {
    return mMode;
  }

  /**
   * Get the directory.
   * @return the directory.
   */
  @Override
  public File directory() {
    return mReaderParams.directory();
  }

  /**
   * Get the first sequence to be processed.
   * @return the first sequence to be processed.
   */
  @Override
  public HashingRegion region() {
    return mRegion;
  }

  /**
   * Get the first sequence to be processed.
   * @return the first sequence to be processed.
   */
  @Override
  public LongRange readerRestriction() {
    return mReaderRestriction;
  }

  /**
   * Get the total number of sequences.
   * @return the total number of sequences.
   */
  @Override
  public long numberSequences() {
    if (mRegion != HashingRegion.NONE) {
      return mRegion.getEnd() - mRegion.getStart();
    } else {
      return reader().numberSequences();
    }
  }

  /**
   * Get a SequencesReader for this sequence.
   * @return a SequencesReader for this sequence. A single reader is returned on succesive calls.
   */
  @Override
  public SequencesReader reader() {
    return mReaderParams.reader();
  }

  /**
   *  Get the reader parameters that specify the directory.
   * @return the reader parameters that specify the directory.
   */
  @Override
  public ReaderParams readerParams() {
    return mReaderParams;
  }

  /**
   * Get the length of the longest sequence.
   * @return the length of the longest sequence.
   */
  @Override
  public long maxLength() {
    return mReaderParams.maxLength();
  }

  /**
   * If necessary carefully close the reader.
   * @throws IOException if an IO error occurs
   */
  @Override
  public void close() throws IOException {
    mReaderParams.close();
  }

  @Override
  public int hashCode() {
    return Utils.pairHash(Utils.pairHash(mReaderParams.hashCode(), (int) mRegion.getStart()), (int) mRegion.getEnd());
  }

  @Override
  public boolean equals(final Object obj) {
    if (this == obj) {
      return true;
    }
    if (obj == null) {
      return false;
    }
    final MockSequenceParams that = (MockSequenceParams) obj;
    if (!this.mReaderParams.equals(that.mReaderParams)) {
      return false;
    }
    if (!this.mRegion.equals(that.mRegion)) {
      return false;
    }
    return true;
  }


  @Override
  public void toString(final StringBuilder sb) {
    toString("", sb);
  }

  void toString(final String prefix, final StringBuilder sb) {
    sb.append(prefix).append("SequenceParams mode=").append(mode()).append(" region=").append(mRegion).append(" directory=").append(mReaderParams.directory());
  }

  @Override
  public boolean integrity() {
    if (mReaderParams == null) {
      throw new NullPointerException();
    }
    if (mRegion != HashingRegion.NONE) {
      Exam.assertTrue(0 <= mRegion.getStart() && mRegion.getStart() <= mRegion.getEnd());
      Exam.assertTrue("end=" + mRegion.getEnd() + " number=" + reader().numberSequences(), mRegion.getEnd() <= reader().numberSequences());
    }
    return true;
  }

  @Override
  public boolean globalIntegrity() {
    integrity();
    Exam.integrity(mReaderParams);
    return true;
  }

    @Override
    public Sex sex() {
      return Sex.FEMALE;
    }

}

