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
package com.rtg.reader;

import java.io.File;
import java.io.IOException;

import com.rtg.mode.SequenceType;

/**
 * Made for a special purpose to behave like two sequence readers are interleaved.
 * If you want more functionality you'll have to implement it.
 */
public class AlternatingSequencesReader implements SequencesReader {

  private final SequencesReader mFirst;
  private final SequencesReader mSecond;

  private long mSequenceId;
  private boolean mCurrentIsFirst;

  /**
   * Constructs a sequence reader which alternates between two given sequence readers.
   * @param first the first sequence reader
   * @param second the second sequence reader
   */
  public AlternatingSequencesReader(final SequencesReader first, final SequencesReader second) {
    mFirst = first;
    mSecond = second;
    mSequenceId = -1;
  }

  @Override
  public long dataChecksum() {
    throw new UnsupportedOperationException("Not supported yet.");
  }
  @Override
  public long qualityChecksum() {
    throw new UnsupportedOperationException("Not supported yet.");
  }
  @Override
  public long nameChecksum() {
    throw new UnsupportedOperationException("Not supported yet.");
  }

  @Override
  public void close() throws IOException {
    try {
      mFirst.close();
    } finally {
      mSecond.close();
    }
  }

  @Override
  public SequencesReader copy() {
    throw new UnsupportedOperationException("Not supported yet.");
  }

  @Override
  public int currentLength() {
    if (mSequenceId < 0) {
      throw new IllegalStateException();
    }
    if (mCurrentIsFirst) {
      return mFirst.currentLength();
    } else {
      return mSecond.currentLength();
    }
  }

  @Override
  public String currentName() throws IllegalStateException, IOException {
    throw new UnsupportedOperationException("Not supported yet.");
  }

  @Override
  public String currentFullName() throws IllegalStateException, IOException {
    return currentName();
  }

  @Override
  public long currentSequenceId() {
    if (mSequenceId < 0) {
      throw new IllegalStateException();
    }
    return mSequenceId;
  }

  @Override
  public File path() {
    throw new UnsupportedOperationException("Not supported yet.");
  }

  @Override
  public PrereadArm getArm() {
    throw new UnsupportedOperationException("Not supported yet.");
  }

  @Override
  public PrereadType getPrereadType() {
    throw new UnsupportedOperationException("Not supported yet.");
  }

  @Override
  public double globalQualityAverage() {
    throw new UnsupportedOperationException("Not supported yet.");
  }

  @Override
  public SdfId getSdfId() {
    throw new UnsupportedOperationException("Not supported yet.");
  }

  @Override
  public boolean hasHistogram() {
    throw new UnsupportedOperationException("Not supported yet.");
  }

  @Override
  public boolean hasQualityData() {
    throw new UnsupportedOperationException("Not supported yet.");
  }

  @Override
  public boolean hasNames() {
    throw new UnsupportedOperationException("Not supported yet.");
  }

  @Override
  public long[] histogram() {
    throw new UnsupportedOperationException("Not supported yet.");
  }

  @Override
  public long lengthBetween(long start, long end) {
    throw new UnsupportedOperationException("Not supported yet.");
  }

  @Override
  public long longestNBlock() {
    throw new UnsupportedOperationException("Not supported yet.");
  }

  @Override
  public long maxLength() {
    throw new UnsupportedOperationException("Not supported yet.");
  }

  @Override
  public long minLength() {
    throw new UnsupportedOperationException("Not supported yet.");
  }

  @Override
  public long nBlockCount() {
    throw new UnsupportedOperationException("Not supported yet.");
  }

  @Override
  public PrereadNames names() {
    throw new UnsupportedOperationException("Not supported yet.");
  }

  @Override
  public boolean nextSequence() throws IOException {
    if (mSequenceId + 1 < numberSequences()) {
      seek(mSequenceId + 1);
      return true;
    }
    return false;
  }

  @Override
  public long numberSequences() {
    return mFirst.numberSequences() + mSecond.numberSequences();
  }

  @Override
  public long[] posHistogram() {
    throw new UnsupportedOperationException("Not supported yet.");
  }

  @Override
  public double[] positionQualityAverage() {
    throw new UnsupportedOperationException("Not supported yet.");
  }

  @Override
  public int length(long sequenceIndex)  {
    throw new UnsupportedOperationException("Not supported yet.");
  }
  @Override
  public String name(long sequenceIndex)  {
    throw new UnsupportedOperationException("Not supported yet.");
  }
  @Override
  public String fullName(long sequenceIndex)  {
    return name(sequenceIndex);
  }
  @Override
  public int read(long sequenceIndex, byte[] dataOut) throws IllegalArgumentException, IOException {
    throw new UnsupportedOperationException("Not supported yet.");
  }
  @Override
  public int read(long sequenceIndex, byte[] dataOut, int start, int length) throws IllegalArgumentException, IllegalStateException, IOException {
    throw new UnsupportedOperationException("Not supported yet.");
  }

  @Override
  public int readCurrent(byte[] dataOut) throws IllegalArgumentException, IllegalStateException, IOException {
    throw new UnsupportedOperationException("Not supported yet.");
  }
  @Override
  public int readCurrent(byte[] dataOut, int start, int length) throws IllegalArgumentException, IllegalStateException, IOException {
    throw new UnsupportedOperationException("Not supported yet.");
  }

  @Override
  public int readCurrentQuality(byte[] dest) throws IllegalArgumentException, IllegalStateException, IOException {
    throw new UnsupportedOperationException("Not supported yet.");
  }

  @Override
  public int readCurrentQuality(byte[] dest, int start, int length) throws IllegalArgumentException, IllegalStateException, IOException {
    throw new UnsupportedOperationException("Not supported yet.");
  }

  @Override
  public int readQuality(long sequenceIndex, byte[] dest) throws IllegalArgumentException, IOException {
    throw new UnsupportedOperationException("Not supported yet.");
  }

  @Override
  public int readQuality(long sequenceIndex, byte[] dest, int start, int length) throws IllegalArgumentException, IOException {
    throw new UnsupportedOperationException("Not supported yet.");
  }

  @Override
  public long[] residueCounts() {
    throw new UnsupportedOperationException("Not supported yet.");
  }

  @Override
  public long sdfVersion() {
    throw new UnsupportedOperationException("Not supported yet.");
  }

  @Override
  public void seek(final long sequenceId) throws IOException {
    mSequenceId = sequenceId;
    mCurrentIsFirst = mSequenceId % 2 == 0;
    if (mCurrentIsFirst) {
      mFirst.seek(mSequenceId / 2);
    } else {
      mSecond.seek(mSequenceId / 2);
    }
  }

  @Override
  public int[] sequenceLengths(long start, long end) {
    throw new UnsupportedOperationException("Not supported yet.");
  }

  @Override
  public long totalLength() {
    return mFirst.totalLength() + mSecond.totalLength();
  }

  @Override
  public SequenceType type() {
    throw new UnsupportedOperationException("Not supported yet.");
  }

  @Override
  public boolean compressed() {
    throw new UnsupportedOperationException("Not supported yet.");
  }

  @Override
  public void reset() {
    mFirst.reset();
    mSecond.reset();
    mSequenceId = -1;
    mCurrentIsFirst = false;
  }

  @Override
  public String currentNameSuffix() throws IllegalStateException, IOException {
    throw new UnsupportedOperationException("Not supported yet.");
  }

  @Override
  public String nameSuffix(long sequenceIndex) {
    throw new UnsupportedOperationException("Not supported yet.");
  }

  @Override
  public long suffixChecksum() {
    throw new UnsupportedOperationException("Not supported yet.");
  }

  @Override
  public String getReadMe() {
    return null;
  }

}
