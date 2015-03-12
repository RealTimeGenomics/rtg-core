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
package com.rtg.launcher;

import java.io.File;
import java.io.IOException;

import com.rtg.mode.SequenceMode;
import com.rtg.reader.MockSequencesReader;
import com.rtg.reader.SequencesReader;
import com.rtg.util.Utils;
import com.rtg.util.intervals.LongRange;

/**
 */
public class MockReaderParams extends ReaderParams {

  private final SequenceMode mMode;

  private final long mMaxLength;

  private final SequencesReader mReader;

  /**
   * Create a new {@link MockReaderParams}
   * @param length mock length for the mock reader.
   * @param numberSequences number of sequences.
   * @param mode {@link SequenceMode}
   */
  public MockReaderParams(final long length, final long numberSequences, final SequenceMode mode) {
    this(new MockSequencesReader(mode.codeType(), numberSequences, length), mode);
  }

  /**
   * Create a new {@link MockReaderParams}
   * @param reader MockReader
   * @param mode {@link SequenceMode}
   */
  public MockReaderParams(final SequencesReader reader, final SequenceMode mode) {
    mReader = reader;
    long mx;
    try {
      mx = reader.maxLength();
    } catch (final UnsupportedOperationException e) {
      mx = -1;
    }
    mMaxLength = mx;
    mMode = mode;
  }

  /**
   * @see ReaderParams#close()
   */
  @Override
  public void close() {
  }

  /**
   * @see ReaderParams#directory()
   */
  @Override
  public File directory() {
    //the serialize-r checks this
    return new File("temp");
  }

  /**
   * @see ReaderParams#maxLength()
   */
  @Override
  public long maxLength() {
    if (mMaxLength == -1) {
      throw new UnsupportedOperationException();
    }
    return mMaxLength;
  }

  /**
   * @see ReaderParams#mode()
   */
  @Override
  public SequenceMode mode() {
    return mMode;
  }

  /**
   * @see ReaderParams#reader()
   */
  @Override
  public SequencesReader reader() {
    return mReader;
  }

  @Override
  public int[] lengths() throws IOException {
    final int n = (int) mReader.numberSequences();
    final int[] lengths = new int[n];
    for (int i = 0; i < n; i++) {
      lengths[i] = mReader.length(i);
    }
    return lengths;
  }

  void toString(final String prefix, final StringBuilder sb) {
    sb.append(prefix).append("SequenceParams mode=").append(mMode).append(" directory=").append(directory());
  }

  @Override
  public int hashCode() {
    return Utils.pairHash(mMode.hashCode(), directory().hashCode());
  }

  @Override
  public boolean equals(final Object obj) {
    if (this == obj) {
      return true;
    }
    if (obj == null) {
      return false;
    }
    final MockReaderParams that = (MockReaderParams) obj;
    return this.mMode.equals(that.mMode) && this.directory().equals(that.directory());
  }

  @Override
  public LongRange adjustedRegion() {
    return null;
  }
}

