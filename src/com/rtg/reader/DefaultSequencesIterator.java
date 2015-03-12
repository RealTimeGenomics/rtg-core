/*
 * Copyright (c) 2015. Real Time Genomics Limited.
 *
 * Use of this source code is bound by the Real Time Genomics Limited Software Licence Agreement
 * for Academic Non-commercial Research Purposes only.
 *
 * If you did not receive a license accompanying this file, a copy must first be obtained by email
 * from support@realtimegenomics.com.  On downloading, using and/or continuing to use this source
 * code you accept the terms of that license agreement and any amendments to those terms that may
 * be made from time to time by Real Time Genomics Limited.
 */
package com.rtg.reader;

import java.io.IOException;

import com.reeltwo.jumble.annotations.TestClass;

/**
 * Provides iterator-style access to a set of sequences provided by an existing SequencesReader.
 */
@TestClass("com.rtg.reader.DefaultSequencesReaderTest")
public class DefaultSequencesIterator implements SequencesIterator {

  private final SequencesReader mReader;
  private final long mNumberSequences;
  private long mSequenceId = -1;

  /**
   * Create the iterator accessor around an existing reader.
   * @param reader the reader
   */
  public DefaultSequencesIterator(SequencesReader reader) {
    mReader = reader;
    mNumberSequences = reader.numberSequences();
  }

  @Override
  public SequencesReader reader() {
    return mReader;
  }

  @Override
  public void reset() {
    mSequenceId = -1;
  }

  @Override
  public void seek(long sequenceId) throws IOException {
    mSequenceId = sequenceId;
    if (mSequenceId >= mNumberSequences || mSequenceId < 0) {
      throw new IllegalArgumentException("Failed to seek to sequence: " + sequenceId + " numberOfSequences: " + mNumberSequences);
    }
  }

  @Override
  public boolean nextSequence() throws IOException {
    return ++mSequenceId < mNumberSequences;
  }

  private void checkSequenceId() {
    if (mSequenceId >= mNumberSequences || mSequenceId < 0) {
      throw new IllegalStateException("Last call to nextSequence() or seek() failed and left current information unavailable.");
    }
  }

  @Override
  public long currentSequenceId() throws IllegalStateException {
    checkSequenceId();
    return mSequenceId;
  }

  @Override
  public int currentLength() throws IllegalStateException, IOException {
    checkSequenceId();
    return mReader.length(mSequenceId);
  }

  @Override
  public String currentName() throws IllegalStateException, IOException {
    checkSequenceId();
    return mReader.name(mSequenceId);
  }

  @Override
  public String currentFullName() throws IllegalStateException, IOException {
    checkSequenceId();
    return mReader.fullName(mSequenceId);
  }

  @Override
  public String currentNameSuffix() throws IllegalStateException, IOException {
    checkSequenceId();
    return mReader.nameSuffix(mSequenceId);
  }


  @Override
  public int readCurrent(byte[] dataOut) throws IllegalArgumentException, IllegalStateException, IOException {
    checkSequenceId();
    return mReader.read(mSequenceId, dataOut);
  }

  @Override
  public int readCurrent(byte[] dataOut, int start, int length) throws IllegalArgumentException, IOException {
    checkSequenceId();
    return mReader.read(mSequenceId, dataOut, start, length);
  }

  @Override
  public int readCurrentQuality(byte[] dest) throws IllegalArgumentException, IllegalStateException, IOException {
    checkSequenceId();
    return mReader.readQuality(mSequenceId, dest);
  }

  @Override
  public int readCurrentQuality(byte[] dest, int start, int length) throws IllegalArgumentException, IllegalStateException, IOException {
    checkSequenceId();
    return mReader.readQuality(mSequenceId, dest, start, length);
  }

}
